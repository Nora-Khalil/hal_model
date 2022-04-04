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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01883,0.0509553,-4.0285e-05,1.75587e-08,-3.55112e-12,-34563.2,22.0219], Tmin=(100,'K'), Tmax=(1054.46,'K')), NASAPolynomial(coeffs=[6.07771,0.0355584,-1.83825e-05,3.71126e-09,-2.68088e-13,-35419.2,2.22484], Tmin=(1054.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(CsCdCsF1s)"""),
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
    label = '[O]C([O])=CC1(F)C=C1F(14260)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {9,S}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {5,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-168.699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,323,467,575,827,1418,3010,987.5,1337.5,450,1655,350,440,435,1725,323.426,323.426,323.426,323.426],'cm^-1')),
        HinderedRotor(inertia=(0.313243,'amu*angstrom^2'), symmetry=1, barrier=(23.2519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762642,0.0743269,-9.26933e-05,5.84308e-08,-1.46301e-11,-20175.8,25.7368], Tmin=(100,'K'), Tmax=(972.713,'K')), NASAPolynomial(coeffs=[13.7265,0.0210178,-1.04881e-05,2.0912e-09,-1.50326e-13,-22697.9,-36.4481], Tmin=(972.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCCCF) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(Cd-Cs(F)-Cd) + radical(C=COJ) + radical(C=COJ) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

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
    label = 'O=C1[CH]C(F)=C[C](F)O1(14261)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {8,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,D}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  C u0 p0 c0 {5,D} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  C u0 p0 c0 {3,S} {4,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-429.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75189,0.0307569,4.47491e-05,-7.3529e-08,2.65644e-11,-51528,16.5308], Tmin=(100,'K'), Tmax=(1091.13,'K')), NASAPolynomial(coeffs=[15.9384,0.0265922,-1.52949e-05,3.34129e-09,-2.55152e-13,-57471.8,-66.1985], Tmin=(1091.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-429.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclohexane) + radical(C=CCJCO) + radical(Cs_P)"""),
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
    label = 'O=[C]C1[C](F)C=C1F(14262)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u1 p0 c0 {3,D} {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-89.2466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,1855,455,950,180,180,180,1054.87,1054.94,1055,1055.01],'cm^-1')),
        HinderedRotor(inertia=(0.00686835,'amu*angstrom^2'), symmetry=1, barrier=(5.4256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61954,0.0492405,-4.22839e-05,1.79364e-08,-3.04945e-12,-10645.8,23.1221], Tmin=(100,'K'), Tmax=(1393.41,'K')), NASAPolynomial(coeffs=[12.7056,0.017416,-8.02501e-06,1.54542e-09,-1.08646e-13,-13735.3,-34.0398], Tmin=(1393.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.2466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=C1OC2(F)C=C(F)C12(14211)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-508.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68647,0.0491332,-3.19732e-05,9.00936e-09,-9.7827e-13,-61102.2,16.539], Tmin=(100,'K'), Tmax=(2127.57,'K')), NASAPolynomial(coeffs=[19.493,0.015655,-8.3697e-06,1.61317e-09,-1.09168e-13,-68679,-82.8106], Tmin=(2127.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-508.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFO) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'O=C(O)C1=C(F)C=C1F(14206)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {9,S} {11,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {5,D} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-266.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.974606,0.0720833,-9.49097e-05,6.53681e-08,-1.81891e-11,-31967.6,24.3348], Tmin=(100,'K'), Tmax=(871.712,'K')), NASAPolynomial(coeffs=[11.3677,0.0243917,-1.28422e-05,2.60305e-09,-1.88144e-13,-33779.5,-24.3785], Tmin=(871.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(CdCCF) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cd-Cd-Cd-Cd(F))"""),
)

species(
    label = 'O=C(O)C1C(F)=C=C1F(14214)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {8,S} {11,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u0 p0 c0 {6,D} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-248.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2419,0.0672989,-9.95705e-05,8.99222e-08,-3.35304e-11,-29748.2,24.0252], Tmin=(100,'K'), Tmax=(756.593,'K')), NASAPolynomial(coeffs=[5.93988,0.0342574,-1.77985e-05,3.53769e-09,-2.5088e-13,-30224.3,4.22235], Tmin=(756.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CdCddCF) + group(CdCddCF) + group(Cds-OdCsOs) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[O]C(=O)C1C2(F)[CH]C12F(14263)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-157.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32735,0.0688421,-0.000114155,1.15981e-07,-4.64841e-11,-18913.8,24.8698], Tmin=(100,'K'), Tmax=(787.842,'K')), NASAPolynomial(coeffs=[2.46201,0.0419486,-2.27164e-05,4.5598e-09,-3.23735e-13,-18436.8,23.8285], Tmin=(787.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCCF) + group(CsCCCF) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-OdCsOs) + polycyclic(s2_3_3_ane) + radical(CCOJ) + radical(bicyclo[1.1.0]butane-secondary) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1C=C(F)C1[C]1OO1(14264)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {9,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (73.1864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.134,0.0452921,9.56676e-06,-5.88158e-08,2.96297e-11,8921.83,26.3336], Tmin=(100,'K'), Tmax=(933.435,'K')), NASAPolynomial(coeffs=[20.8022,0.00618234,-1.76777e-07,-1.13413e-11,-5.4443e-15,3282.07,-77.7407], Tmin=(933.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.1864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFH) + group(Cs-CsOsOsH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(dioxirane) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cs_P)"""),
)

species(
    label = 'O=C1OC2(F)[CH][C](F)C12(14239)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7  C u1 p0 c0 {2,S} {5,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-235.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70752,0.044759,-2.23128e-05,1.64513e-09,8.80797e-13,-28269.3,22.1298], Tmin=(100,'K'), Tmax=(1476.96,'K')), NASAPolynomial(coeffs=[14.5537,0.0220154,-1.14497e-05,2.26452e-09,-1.58858e-13,-33378,-49.3039], Tmin=(1476.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-235.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFO) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-OdCsOs) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(CCJCO)"""),
)

species(
    label = 'O=C1OC2[C](F)C1[C]2F(14222)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {11,S}
7  C u1 p0 c0 {1,S} {5,S} {6,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-277.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0711,0.0406737,-1.88265e-05,1.43534e-09,5.32463e-13,-33303.1,21.4326], Tmin=(100,'K'), Tmax=(1664.51,'K')), NASAPolynomial(coeffs=[14.5192,0.0213618,-1.09777e-05,2.11846e-09,-1.4489e-13,-38915.8,-49.3772], Tmin=(1664.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cds-OdCsOs) + polycyclic(s3_4_5_ane) + radical(CsCsCsF1s) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(=O)C=C(F)C=[C]F(14265)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {9,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u0 p0 c0 {5,D} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {3,S} {4,D} {6,S}
9  C u1 p0 c0 {2,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-141.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2995,3025,975,1000,1300,1375,400,500,1630,1680,167,640,1190,180,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.33013,'amu*angstrom^2'), symmetry=1, barrier=(7.59033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329466,'amu*angstrom^2'), symmetry=1, barrier=(7.57506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24234,0.0740466,-0.000148331,1.50881e-07,-5.7186e-11,-16906.2,24.2187], Tmin=(100,'K'), Tmax=(837.353,'K')), NASAPolynomial(coeffs=[3.76248,0.0315892,-1.77836e-05,3.56025e-09,-2.49337e-13,-16261.8,18.8757], Tmin=(837.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CdCCF) + missing(Cd-COCdH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFH) + radical(CCOJ) + radical(Cdj(Cd-CdH)(F1s))"""),
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
    label = 'O=C(O)[C]1[C](F)C=C1F(14267)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {9,S} {11,S}
4  O u0 p2 c0 {9,D}
5  C u1 p0 c0 {6,S} {7,S} {9,S}
6  C u1 p0 c0 {1,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {10,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,346,659,817,1284,271,519,563,612,1379,2950,1000,361.806,964.567,964.567,964.567,964.567,964.567,964.567,964.567,964.567,964.567,2329.39],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33218,0.0541601,-4.1717e-05,1.53714e-08,-2.26361e-12,-43324.6,24.5058], Tmin=(100,'K'), Tmax=(1585,'K')), NASAPolynomial(coeffs=[14.8177,0.0201273,-9.50918e-06,1.82445e-09,-1.26866e-13,-47599.4,-46.7651], Tmin=(1585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-361.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCJ(C)CO) + radical(CsCdCsF1s)"""),
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
    label = '[O]C(=O)C1=C(F)C=C1F(14268)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,D} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {5,D} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {10,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-40.9547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([217,343,449,587,660,812,790,914,798,948,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03313,0.0732738,-0.000122853,1.13959e-07,-4.20399e-11,-4826.65,24.2697], Tmin=(100,'K'), Tmax=(778.936,'K')), NASAPolynomial(coeffs=[7.14962,0.0299601,-1.65192e-05,3.33123e-09,-2.36859e-13,-5418.38,-1.39262], Tmin=(778.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.9547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(CdCCF) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cd-Cd-Cd-Cd(F)) + radical(CCOJ)"""),
)

species(
    label = '[O]C(=O)C1C(F)=C=C1F(14269)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u0 p0 c0 {6,D} {7,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-22.4083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,104,186,245,407,330,466,761,907,1218,1388,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52896,0.0654971,-0.000115424,1.1992e-07,-4.78113e-11,-2617.02,23.1581], Tmin=(100,'K'), Tmax=(811.537,'K')), NASAPolynomial(coeffs=[1.49379,0.0402361,-2.17214e-05,4.32559e-09,-3.04648e-13,-1773.76,28.4807], Tmin=(811.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.4083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CdCddCF) + group(CdCddCF) + group(Cds-OdCsOs) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(CCOJ)"""),
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
    label = '[O]C(=O)C1=C(F)C=[C]1(14270)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {7,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {4,D} {6,S}
6 C u0 p0 c0 {5,S} {8,D} {9,S}
7 C u0 p0 c0 {2,S} {3,D} {4,S}
8 C u1 p0 c0 {4,S} {6,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (373.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2950,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36542,0.0656676,-0.000114205,1.05883e-07,-3.812e-11,45017.3,23.6755], Tmin=(100,'K'), Tmax=(819.861,'K')), NASAPolynomial(coeffs=[6.85773,0.0246632,-1.31901e-05,2.60518e-09,-1.8214e-13,44594.2,1.18129], Tmin=(819.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cd-Cd-Cd-Cd(F)) + radical(CCOJ) + radical(cyclobutadiene-C1)"""),
)

species(
    label = '[O]C(=O)C1[C]=C=C1F(14271)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5 C u0 p0 c0 {1,S} {4,S} {8,D}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 C u1 p0 c0 {4,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (406.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,145,326,398,834,1303,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68261,0.0268827,-1.12461e-05,8.38551e-10,1.18045e-13,48809.9,12.4213], Tmin=(100,'K'), Tmax=(2721.77,'K')), NASAPolynomial(coeffs=[43.9971,-0.0145777,2.61046e-06,-3.5278e-10,2.51507e-14,21364.8,-227.716], Tmin=(2721.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CdCddCF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C(=O)[C]1C(F)=CC1F(14272)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u0 p0 c0 {3,S} {4,D} {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-273.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2950,1000,271,519,563,612,1379,600,879.515,879.515,879.515,879.515,879.515,879.515,879.515,879.515,879.515,2418.22],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61196,0.0428532,-2.38916e-05,5.3938e-09,-4.52357e-13,-32832.1,20.5332], Tmin=(100,'K'), Tmax=(2834.85,'K')), NASAPolynomial(coeffs=[27.8935,0.00718255,-5.01818e-06,9.55582e-10,-6.09786e-14,-47166.7,-127.78], Tmin=(2834.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[O]C(=O)C1C(F)=[C]C1F(14273)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-161.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,164,312,561,654,898,1207,1299,3167,246,474,533,1155,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59313,0.0586369,-7.30403e-05,6.34299e-08,-2.48759e-11,-19290.6,23.2647], Tmin=(100,'K'), Tmax=(668.995,'K')), NASAPolynomial(coeffs=[4.61331,0.0377314,-1.97823e-05,3.99499e-09,-2.87824e-13,-19630.9,10.3844], Tmin=(668.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCOJ) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
)

species(
    label = 'O=C(O)C1[C](F)[C]=C1F(14274)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {8,S} {11,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-248.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,1000,346,659,817,1284,246,474,533,1155,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35592,0.0583542,-5.22269e-05,2.34409e-08,-4.28676e-12,-29831.9,24.1255], Tmin=(100,'K'), Tmax=(1283.74,'K')), NASAPolynomial(coeffs=[12.5617,0.0234382,-1.14286e-05,2.25364e-09,-1.60654e-13,-32709,-32.7346], Tmin=(1283.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
)

species(
    label = 'O=C(OF)C1[C]=C[C]1F(14275)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {4,D} {5,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (81.8704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,3150,900,1100,346,659,817,1284,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5076,0.0582627,-5.17892e-05,2.34677e-08,-4.43616e-12,9933.44,24.2222], Tmin=(100,'K'), Tmax=(1220.46,'K')), NASAPolynomial(coeffs=[10.9048,0.0274633,-1.39349e-05,2.78978e-09,-2.00433e-13,7639.68,-22.9861], Tmin=(1220.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.8704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[O]C(=O)C1[C]=CC1(F)F(14276)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {5,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-192.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,274,345,380,539,705,1166,1213,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53247,0.0601617,-7.87096e-05,7.05795e-08,-2.78351e-11,-23070.5,23.4893], Tmin=(100,'K'), Tmax=(697.521,'K')), NASAPolynomial(coeffs=[4.57074,0.0377448,-1.9764e-05,3.97768e-09,-2.8555e-13,-23372.9,10.7966], Tmin=(697.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(CCOJ) + radical(cyclobutene-vinyl)"""),
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
    E0 = (-178.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (110.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-23.4092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-80.2649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (263.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-170.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-100.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-137.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (35.2749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (182.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (63.7824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-37.4107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (106.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (68.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (8.41634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (170.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (280.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (298.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-103.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (227.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (315.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (426.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (19.7083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (131.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (12.3631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (253.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (99.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['CO2(14)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C([O])=CC1(F)C=C1F(14260)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.66547e+10,'s^-1'), n=0.732, Ea=(169.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCsCJ;CsJ-OneDeH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-OneDeH;C]
Euclidian distance = 1.0
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
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['O=C1[CH]C(F)=C[C](F)O1(14261)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.75746e+12,'s^-1'), n=0.324012, Ea=(98.1226,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=[C]C1[C](F)C=C1F(14262)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['O=C1OC2(F)C=C(F)C12(14211)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.48e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 4.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['O=C(O)C1=C(F)C=C1F(14206)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['O=C(O)C1C(F)=C=C1F(14214)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(40.8634,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['[O]C(=O)C1C2(F)[CH]C12F(14263)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['F[C]1C=C(F)C1[C]1OO1(14264)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(361.103,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 359.8 to 361.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['O=C1OC2(F)[CH][C](F)C12(14239)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.012e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['O=C1OC2[C](F)C1[C]2F(14222)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.49279e+11,'s^-1'), n=0.14706, Ea=(140.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra;radadd_intra] for rate rule [Rn2c4_beta;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(=O)C=C(F)C=[C]F(14265)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingle]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    products = ['[O]C1([O])C2C(F)=CC21F(14266)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.27712e+09,'s^-1'), n=0.755479, Ea=(247.37,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C(O)[C]1[C](F)C=C1F(14267)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(259.935,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=O(517)', 'FC1=CC(F)=C1(307)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(306.062,'m^3/(mol*s)'), n=1.16366, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[O]C(=O)C1=C(F)C=C1F(14268)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[O]C(=O)C1C(F)=C=C1F(14269)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.10362,'m^3/(mol*s)'), n=2.13572, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R_Ext-5R!H-R_Ext-2CS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R_Ext-5R!H-R_Ext-2CS-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CO2(14)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.33888e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(103.074,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][C]=O(517)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C(=O)C1=C(F)C=[C]1(14270)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(113.709,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[O]C(=O)C1[C]=C=C1F(14271)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(191.966,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(=O)[C]1C(F)=CC1F(14272)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.64354e+09,'s^-1'), n=1.15628, Ea=(183.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_noH;Cs_H_out_noH] + [R2H_S;C_rad_out_TwoDe;Cs_H_out] for rate rule [R2H_S_cy4;C_rad_out_TwoDe;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)C1C(F)=[C]C1F(14273)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C(O)C1[C](F)[C]=C1F(14274)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C(OF)C1[C]=C[C]1F(14275)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(62.1856,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(=O)C1[C]=CC1(F)F(14276)'],
    products = ['[O]C(=O)C1[C](F)C=C1F(14208)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(182.661,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #4041',
    isomers = [
        '[O]C(=O)C1[C](F)C=C1F(14208)',
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
    label = 'PDepNetwork #4041',
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

