species(
    label = '[CH]=C(F)C([CH]F)=C(F)F(3021)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {8,D}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u1 p0 c0 {2,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {5,D}
9  C u1 p0 c0 {6,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-303.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,250,446,589,854,899,234,589,736,816,1240,3237,182,240,577,636,1210,1413,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.163522,'amu*angstrom^2'), symmetry=1, barrier=(13.5131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07233,'amu*angstrom^2'), symmetry=1, barrier=(47.647,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507821,0.0815657,-0.000109325,7.40995e-08,-1.99383e-11,-36378.3,25.3329], Tmin=(100,'K'), Tmax=(907.429,'K')), NASAPolynomial(coeffs=[13.9112,0.0224817,-1.16554e-05,2.34252e-09,-1.68532e-13,-38810.7,-38.0282], Tmin=(907.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(CdCFF) + group(Cds-CdsHH) + radical(CsCdF1sH) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[CH]=C(F)C(F)[C]=C(F)F(3022)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {9,D}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u1 p0 c0 {5,S} {7,D}
9  C u1 p0 c0 {6,D} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-208.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,246,474,533,1155,562,600,623,1070,1265,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.391366,'amu*angstrom^2'), symmetry=1, barrier=(8.99827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11365,'amu*angstrom^2'), symmetry=1, barrier=(25.605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3183.59,'J/mol'), sigma=(5.0267,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=497.27 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462081,0.0870705,-0.000149319,1.34707e-07,-4.72222e-11,-24979.9,27.9096], Tmin=(100,'K'), Tmax=(826.507,'K')), NASAPolynomial(coeffs=[8.99023,0.0295731,-1.55243e-05,3.03665e-09,-2.10961e-13,-25835.5,-8.25639], Tmin=(826.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=C=C([CH]F)C(F)(F)F(3604)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u1 p0 c0 {4,S} {6,S} {10,S}
8  C u0 p0 c0 {6,D} {9,D}
9  C u1 p0 c0 {8,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-415.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518696,0.080896,-0.000113808,8.28542e-08,-2.38962e-11,-49859.4,25.7344], Tmin=(100,'K'), Tmax=(850.195,'K')), NASAPolynomial(coeffs=[12.9331,0.0224941,-1.0779e-05,2.07337e-09,-1.44795e-13,-51970.5,-32.1443], Tmin=(850.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-415.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCdF1sH) + radical(C=C=CJ)"""),
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
    label = '[CH]C(F)=C=C(F)F(3605)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 C u2 p0 c0 {4,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-48.3714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([145,326,398,834,1303,94,120,354,641,825,1294,540,610,2055,357.174,357.267,358.75],'cm^-1')),
        HinderedRotor(inertia=(0.00130677,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86031,0.0500541,-5.39129e-05,3.23242e-08,-8.07825e-12,-5743.2,20.2738], Tmin=(100,'K'), Tmax=(953.576,'K')), NASAPolynomial(coeffs=[8.41796,0.0225462,-1.0642e-05,2.07229e-09,-1.47002e-13,-6993.83,-11.0512], Tmin=(953.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.3714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = 'FC1=CC(F)C1=C(F)F(3431)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {9,D}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u0 p0 c0 {3,S} {4,S} {6,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-570.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06062,0.0608518,-5.5166e-05,2.40173e-08,-4.14186e-12,-68533.4,19.1012], Tmin=(100,'K'), Tmax=(1384.16,'K')), NASAPolynomial(coeffs=[15.7912,0.0182828,-9.03434e-06,1.79846e-09,-1.28802e-13,-72611.3,-56.7542], Tmin=(1384.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-570.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C(F)[C]1C(F)C1(F)F(3606)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u1 p0 c0 {8,D} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-279.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29517,0.0611153,-6.0564e-05,3.11432e-08,-6.54663e-12,-33532,22.9311], Tmin=(100,'K'), Tmax=(1129.98,'K')), NASAPolynomial(coeffs=[11.7015,0.0242777,-1.16635e-05,2.29266e-09,-1.63615e-13,-35883.8,-28.5451], Tmin=(1129.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(Allyl_T) + radical(Cdj(Cd-CsF1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH]C1=C(F)[CH]C1(F)F(3607)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u1 p0 c0 {5,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {6,D} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-399.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.833885,0.0672132,-6.87732e-05,3.44324e-08,-6.80875e-12,-47902.6,25.6657], Tmin=(100,'K'), Tmax=(1222.61,'K')), NASAPolynomial(coeffs=[15.7661,0.0183595,-8.83513e-06,1.74911e-09,-1.25613e-13,-51553.8,-49.3746], Tmin=(1222.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(Csj(Cs-F1sF1sCd)(Cd-CdF1s)(H)_ring) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1([C](F)F)C=C1F(3608)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u1 p0 c0 {2,S} {5,S} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-193.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943986,0.0699103,-7.70914e-05,4.21832e-08,-9.25533e-12,-23144.7,26.8238], Tmin=(100,'K'), Tmax=(1095.22,'K')), NASAPolynomial(coeffs=[14.0075,0.022199,-1.1746e-05,2.40676e-09,-1.75701e-13,-26006.2,-37.388], Tmin=(1095.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(CsCsFHH) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd-Cd-Cs(C-FF)) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F(37)',
    structure = adjacencyList("""multiplicity 2
1 F u1 p3 c0
"""),
    E0 = (72.8916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.9984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CC([CH]F)=C(F)F(3609)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u1 p0 c0 {1,S} {4,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {4,D}
7  C u0 p0 c0 {4,S} {8,T}
8  C u0 p0 c0 {7,T} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-197.607,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.44013,'amu*angstrom^2'), symmetry=1, barrier=(38.9766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.699861,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27157,0.0630588,-7.75357e-05,4.96557e-08,-1.27676e-11,-23671,21.5802], Tmin=(100,'K'), Tmax=(944.228,'K')), NASAPolynomial(coeffs=[11.3091,0.0205369,-9.98483e-06,1.96138e-09,-1.39657e-13,-25566.5,-26.2688], Tmin=(944.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCtCs) + group(CdCFF) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CsCdF1sH)"""),
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
    label = 'C#C[C]([C](F)F)C(F)F(3610)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u1 p0 c0 {5,S} {7,S} {8,S}
7  C u1 p0 c0 {3,S} {4,S} {6,S}
8  C u0 p0 c0 {6,S} {9,T}
9  C u0 p0 c0 {8,T} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-366.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0468497,0.0952426,-0.000162539,1.37856e-07,-4.43623e-11,-43949.6,28.5687], Tmin=(100,'K'), Tmax=(898.956,'K')), NASAPolynomial(coeffs=[12.1948,0.0225598,-1.0175e-05,1.81076e-09,-1.17003e-13,-45381,-24.5573], Tmin=(898.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Tert_Propargyl) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'F[CH]C([C]=CF)=C(F)F(3611)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u1 p0 c0 {1,S} {5,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-364.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335557,0.0841789,-0.000118931,8.5599e-08,-2.42051e-11,-43689.3,25.5214], Tmin=(100,'K'), Tmax=(869.875,'K')), NASAPolynomial(coeffs=[14.1229,0.0207796,-9.60611e-06,1.8129e-09,-1.25134e-13,-46087.9,-39.0724], Tmin=(869.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + group(CdCFH) + radical(CsCdF1sH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(F)C(=[C]F)C(F)F(3612)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u1 p0 c0 {4,S} {6,D}
9  C u1 p0 c0 {7,D} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-204.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,250,446,589,854,899,167,640,1190,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.290634,'amu*angstrom^2'), symmetry=1, barrier=(6.68224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867037,'amu*angstrom^2'), symmetry=1, barrier=(19.9349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331334,0.0871246,-0.000134183,1.06782e-07,-3.37453e-11,-24505.6,26.9167], Tmin=(100,'K'), Tmax=(776.422,'K')), NASAPolynomial(coeffs=[12.5271,0.0242913,-1.27877e-05,2.54347e-09,-1.80157e-13,-26399.3,-28.8342], Tmin=(776.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(CdCFH) + group(Cds-CdsHH) + radical(CdCdF1s) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = 'F[C]C(=CF)C(F)=CF(3613)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {5,D} {10,S}
8  C u0 p0 c0 {3,S} {6,D} {11,S}
9  C u2 p0 c0 {4,S} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-302.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,280,518,736,852,873,151,237,609,755,844,966,1147,1245,1323,1443,3181,3261,253.326,253.91,254.279,2025.46],'cm^-1')),
        HinderedRotor(inertia=(0.146653,'amu*angstrom^2'), symmetry=1, barrier=(6.70834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667914,'amu*angstrom^2'), symmetry=1, barrier=(30.8014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516519,0.0833482,-0.000124568,9.79298e-08,-3.0894e-11,-36244.9,24.5722], Tmin=(100,'K'), Tmax=(774.658,'K')), NASAPolynomial(coeffs=[11.5669,0.026286,-1.40711e-05,2.83251e-09,-2.02549e-13,-37956.9,-25.9177], Tmin=(774.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)[C]=CF(3023)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {3,S} {5,S} {9,D}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  C u1 p0 c0 {6,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-228.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,246,474,533,1155,615,860,1140,1343,3152,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27723,'amu*angstrom^2'), symmetry=1, barrier=(29.366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28015,'amu*angstrom^2'), symmetry=1, barrier=(29.4333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3169.09,'J/mol'), sigma=(5.31747,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.00 K, Pc=47.83 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251352,0.0921736,-0.000161104,1.44452e-07,-4.99525e-11,-27391,27.4347], Tmin=(100,'K'), Tmax=(837.248,'K')), NASAPolynomial(coeffs=[9.99882,0.0280464,-1.47577e-05,2.87456e-09,-1.98634e-13,-28407.8,-14.1846], Tmin=(837.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    label = '[CH]C(F)=C=CF(536)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {5,D} {6,S}
4 C u0 p0 c0 {2,S} {5,D} {7,S}
5 C u0 p0 c0 {3,D} {4,D}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (149.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([145,326,398,834,1303,113,247,382,1207,3490,540,610,2055,526.456,526.468,526.48,526.56],'cm^-1')),
        HinderedRotor(inertia=(0.000608275,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32587,0.0400923,-3.38632e-05,1.64559e-08,-3.54973e-12,18091.8,17.0495], Tmin=(100,'K'), Tmax=(1047.8,'K')), NASAPolynomial(coeffs=[6.41226,0.0244923,-1.15305e-05,2.24654e-09,-1.59417e-13,17235.5,-2.85571], Tmin=(1047.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = 'FC=C1C(F)=CC1(F)F(3199)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,D}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {3,S} {6,S} {7,D}
9  C u0 p0 c0 {4,S} {6,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-597.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.923633,0.0641232,-6.16772e-05,2.8543e-08,-5.21154e-12,-71706.1,20.059], Tmin=(100,'K'), Tmax=(1317.09,'K')), NASAPolynomial(coeffs=[16.2267,0.017648,-8.74772e-06,1.75187e-09,-1.26255e-13,-75737.2,-57.9841], Tmin=(1317.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFH) + ring(Cyclobutene)"""),
)

species(
    label = 'F[C](F)[C]1C(F)=CC1F(3614)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-380.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36109,0.0576966,-5.23227e-05,2.40397e-08,-4.49168e-12,-45628.8,22.8464], Tmin=(100,'K'), Tmax=(1263.04,'K')), NASAPolynomial(coeffs=[12.4126,0.022697,-1.07569e-05,2.1002e-09,-1.49092e-13,-48420.5,-33.0516], Tmin=(1263.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCCFH) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Allyl_T) + radical(CsCsF1sF1s)"""),
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
    label = 'C=C(F)C([C]F)=C(F)F(3615)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  C u2 p0 c0 {4,S} {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-336.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742941,0.0778017,-0.000106328,7.6118e-08,-2.19839e-11,-40360.8,24.0022], Tmin=(100,'K'), Tmax=(841.611,'K')), NASAPolynomial(coeffs=[11.672,0.0258586,-1.3751e-05,2.78582e-09,-2.00928e-13,-42200.4,-26.84], Tmin=(841.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(CdCFF) + group(Cds-CdsHH) + radical(CsCCl_triplet)"""),
)

species(
    label = '[CH]=C(F)C(=[CH])C(F)(F)F(3616)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u0 p0 c0 {4,S} {6,S} {9,D}
8  C u1 p0 c0 {6,D} {10,S}
9  C u1 p0 c0 {7,D} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-265.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,350,440,435,1725,250,446,589,854,899,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(0.489674,'amu*angstrom^2'), symmetry=1, barrier=(11.2586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.963027,'amu*angstrom^2'), symmetry=1, barrier=(22.1419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51524,0.0778396,-9.69273e-05,5.92566e-08,-1.41873e-11,-31748.8,25.6276], Tmin=(100,'K'), Tmax=(1022.82,'K')), NASAPolynomial(coeffs=[15.9065,0.0176474,-8.65238e-06,1.71892e-09,-1.23643e-13,-34897.3,-48.9734], Tmin=(1022.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-265.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = '[CH]C(C(F)=CF)=C(F)F(3617)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {3,S} {4,S} {5,D}
8  C u0 p0 c0 {2,S} {6,D} {10,S}
9  C u2 p0 c0 {5,S} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-327.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.62748,0.0784079,-9.18521e-05,5.74585e-08,-1.46664e-11,-39303.6,25.0277], Tmin=(100,'K'), Tmax=(943.971,'K')), NASAPolynomial(coeffs=[12.1121,0.0297423,-1.45198e-05,2.84285e-09,-2.01873e-13,-41471.8,-29.7164], Tmin=(943.971,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFF) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(AllylJ2_triplet)"""),
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
    E0 = (-118.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (36.9731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (103.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (317.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-144.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (43.1851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-14.2446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-42.1549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (136.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (92.0007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (45.8433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (300.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (241.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (45.5464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1.04371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (127.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (90.3627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (48.2884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (334.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-144.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-14.2446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (97.4112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (251.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (95.3429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (63.9211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['C2HF(58)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(34.1072,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 34.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(F)C(F)[C]=C(F)F(3022)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['[CH]=C=C([CH]F)C(F)(F)F(3604)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.72966e+11,'s^-1'), n=0.63878, Ea=(255.44,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]F(804)', '[CH]C(F)=C=C(F)F(3605)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['FC1=CC(F)C1=C(F)F(3431)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_1H;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['[CH]=C(F)[C]1C(F)C1(F)F(3606)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['F[CH]C1=C(F)[CH]C1(F)F(3607)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['F[CH]C1([C](F)F)C=C1F(3608)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(110.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 108.7 to 110.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]F(489)', 'FC=C=C(F)F(1325)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.09106e-06,'m^3/(mol*s)'), n=3.54229, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13159922714560035, var=1.089943264098092, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-2C-R_N-Sp-5R!H-1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-2C-R_N-Sp-5R!H-1R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'C#CC([CH]F)=C(F)F(3609)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(65.5384,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C2HF(58)', 'F[CH][C]=C(F)F(2842)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]F(489)', 'F[CH][C]=C(F)F(2842)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -30.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CHF(40)', '[CH]C(F)=C=C(F)F(3605)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['C#C[C]([C](F)F)C(F)F(3610)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(197.846,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['F[CH]C([C]=CF)=C(F)F(3611)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(153.344,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(F)C(=[C]F)C(F)F(3612)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(181.22,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C]C(=CF)C(F)=CF(3613)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(241.533,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(F)C(F)(F)[C]=CF(3023)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C]F(138)', '[CH]C(F)=C=CF(536)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['FC=C1C(F)=CC1(F)F(3199)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_noH;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['F[C](F)[C]1C(F)=CC1F(3614)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CF2(43)', '[CH]C(F)=C=CF(536)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['C=C(F)C([C]F)=C(F)F(3615)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_single] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_single]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(F)C(=[CH])C(F)(F)F(3616)'],
    products = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(209.176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(F)C([CH]F)=C(F)F(3021)'],
    products = ['[CH]C(C(F)=CF)=C(F)F(3617)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(216.221,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

network(
    label = 'PDepNetwork #1176',
    isomers = [
        '[CH]=C(F)C([CH]F)=C(F)F(3021)',
    ],
    reactants = [
        ('C2HF(58)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1176',
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

