species(
    label = 'F[C]=CC1[C](F)C2OOC12F(1020)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {13,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-45.7609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612869,0.078084,-6.68264e-05,2.75182e-08,-4.62682e-12,-5385.44,26.1816], Tmin=(100,'K'), Tmax=(1369.18,'K')), NASAPolynomial(coeffs=[15.9585,0.0332523,-1.7711e-05,3.6034e-09,-2.60168e-13,-9587.61,-52.6738], Tmin=(1369.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.7609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'FC1=CC2(F)OOC12(939)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-248.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,323,467,575,827,1418,556.295,556.354,556.401,556.506,556.667,556.693,556.734,556.988,556.993],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3642.96,'J/mol'), sigma=(5.89185,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.02 K, Pc=40.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37505,0.0476133,-2.48125e-05,-2.36513e-09,3.49961e-12,-29797.1,16.5191], Tmin=(100,'K'), Tmax=(1189.18,'K')), NASAPolynomial(coeffs=[15.9615,0.016946,-9.33452e-06,1.96667e-09,-1.45909e-13,-34567,-61.8485], Tmin=(1189.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'F[C]=CC1[C](F)C(F)C12OO2(1153)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-38.3553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,259,529,569,1128,1321,1390,3140,212,367,445,1450,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.670341,0.0891496,-7.79888e-05,2.40762e-08,1.79326e-13,-4433.02,32.3668], Tmin=(100,'K'), Tmax=(1023.06,'K')), NASAPolynomial(coeffs=[24.7092,0.0150153,-6.08861e-06,1.20007e-09,-8.96942e-14,-10939.3,-97.0719], Tmin=(1023.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.3553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s1_3_4_ane) + radical(CsCsCsF1s) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]C=CC1(F)C2OOC21F(1154)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {13,S}
10 C u0 p0 c0 {9,D} {11,S} {14,S}
11 C u2 p0 c0 {3,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-112.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605574,0.0838361,-8.27757e-05,4.15831e-08,-8.77691e-12,-13460.7,24.3239], Tmin=(100,'K'), Tmax=(1096.67,'K')), NASAPolynomial(coeffs=[12.8201,0.0392861,-2.1843e-05,4.54307e-09,-3.33412e-13,-16139.9,-35.7316], Tmin=(1096.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_4_ane) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=CC1(F)[CH]C2(F)OOC12(1019)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-54.1072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.608722,0.0840208,-8.3122e-05,4.24426e-08,-9.1657e-12,-6393.32,24.8855], Tmin=(100,'K'), Tmax=(1070.61,'K')), NASAPolynomial(coeffs=[12.2024,0.0407043,-2.24321e-05,4.65092e-09,-3.4086e-13,-8875.78,-31.8384], Tmin=(1070.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.1072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=CC1C2(F)[CH]OOC12F(1155)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {7,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-98.5193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400847,0.0794999,-6.92927e-05,2.93171e-08,-5.02014e-12,-11720.3,27.4487], Tmin=(100,'K'), Tmax=(1363.83,'K')), NASAPolynomial(coeffs=[16.9812,0.0308712,-1.58087e-05,3.17315e-09,-2.27766e-13,-16242.9,-57.6867], Tmin=(1363.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.5193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_3_5_ane) + radical(CCsJOOC) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1=CC2C3(F)OOC3C12F(1156)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
9  C u0 p0 c0 {5,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-346.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784228,0.073422,-5.59967e-05,1.93905e-08,-2.68517e-12,-41617.1,18.1302], Tmin=(100,'K'), Tmax=(1638.98,'K')), NASAPolynomial(coeffs=[18.4349,0.0303447,-1.65722e-05,3.35429e-09,-2.39107e-13,-47402.9,-75.7449], Tmin=(1638.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + polycyclic(s2_4_4_ene_1) + polycyclic(s2_4_4_ane) - ring(Cyclobutane)"""),
)

species(
    label = 'FC=CC1=C(F)C2OOC12F(1157)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,D} {10,S}
9  C u0 p0 c0 {2,S} {7,S} {8,D}
10 C u0 p0 c0 {8,S} {11,D} {13,S}
11 C u0 p0 c0 {3,S} {10,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-390.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.590617,0.0884628,-8.16863e-05,3.50135e-08,-5.81617e-12,-46752.9,25.8627], Tmin=(100,'K'), Tmax=(1460.03,'K')), NASAPolynomial(coeffs=[25.2469,0.0176763,-8.96196e-06,1.8067e-09,-1.30176e-13,-54297.6,-108.567], Tmin=(1460.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-390.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'FC=CC1C(F)=C2OOC21F(1021)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u0 p0 c0 {3,S} {10,D} {14,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-304.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932866,0.080613,-7.0743e-05,-1.27895e-08,5.06653e-11,-36562.7,21.1884], Tmin=(100,'K'), Tmax=(479.957,'K')), NASAPolynomial(coeffs=[7.04845,0.0480323,-2.63838e-05,5.41433e-09,-3.92964e-13,-37361.5,-6.03282], Tmin=(479.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFO) + group(Cds-CdsCsOs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane)"""),
)

species(
    label = 'F[C]=CC1[C](F)OOC=C1F(1158)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {10,D}
8  C u1 p0 c0 {2,S} {4,S} {6,S}
9  C u0 p0 c0 {6,S} {11,D} {13,S}
10 C u0 p0 c0 {5,S} {7,D} {14,S}
11 C u1 p0 c0 {3,S} {9,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-109.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,323,467,575,827,1418,395,473,707,1436,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.179856,0.0979154,-0.000139899,1.10464e-07,-3.51052e-11,-13020,28.3355], Tmin=(100,'K'), Tmax=(794.714,'K')), NASAPolynomial(coeffs=[11.9962,0.0344349,-1.59385e-05,3.00039e-09,-2.06003e-13,-14886,-27.1732], Tmin=(794.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFH) + ring(34dihydro12dioxin) + radical(CsCsF1sO2s) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C]=CC=C(F)C1OO[C]1F(1159)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {7,D} {10,S} {13,S}
10 C u0 p0 c0 {9,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-69.5155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,323,467,575,827,1418,395,473,707,1436,2995,3025,975,1000,1300,1375,400,500,1630,1680,167,640,1190,180,180,180,180,893.153,1331.55,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.199796,'amu*angstrom^2'), symmetry=1, barrier=(4.59371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199796,'amu*angstrom^2'), symmetry=1, barrier=(4.59371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0950068,0.0944258,-0.000118876,7.68969e-08,-1.98291e-11,-8216.96,34.6484], Tmin=(100,'K'), Tmax=(944.456,'K')), NASAPolynomial(coeffs=[15.5287,0.0282565,-1.37862e-05,2.7179e-09,-1.93977e-13,-11168.2,-39.8345], Tmin=(944.456,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.5155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = '[O]OC1(F)C=C(F)C1C=[C]F(1160)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {13,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-136.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,1380,1390,370,380,2900,435,323,467,575,827,1418,3010,987.5,1337.5,450,1655,167,640,1190,272.848,272.849,272.849,272.849,272.849,272.85,840.866,3544.43],'cm^-1')),
        HinderedRotor(inertia=(0.813421,'amu*angstrom^2'), symmetry=1, barrier=(42.972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813415,'amu*angstrom^2'), symmetry=1, barrier=(42.972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.241122,0.100193,-0.000137017,9.97505e-08,-2.92503e-11,-16233.5,29.9171], Tmin=(100,'K'), Tmax=(830.973,'K')), NASAPolynomial(coeffs=[13.5883,0.0336165,-1.68281e-05,3.31735e-09,-2.35402e-13,-18531.6,-34.2401], Tmin=(830.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFO) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=[C]F(252)',
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
    label = 'F[C]=CC1=C(F)C2OOC12F(1161)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,D} {10,S}
9  C u0 p0 c0 {2,S} {7,S} {8,D}
10 C u0 p0 c0 {8,S} {11,D} {13,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-134.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,323,467,575,827,1418,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0508878,0.0889083,-9.23953e-05,4.58518e-08,-9.0067e-12,-16027.7,23.8379], Tmin=(100,'K'), Tmax=(1223.64,'K')), NASAPolynomial(coeffs=[19.6873,0.0243858,-1.33006e-05,2.75944e-09,-2.02599e-13,-20858.2,-75.3711], Tmin=(1223.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + polycyclic(s2_4_4_ene_1) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[C]=CC1C(F)=C2OOC21F(1162)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-54.6891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,323,467,575,827,1418,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0822044,0.0964824,-0.000154779,1.44504e-07,-5.5051e-11,-6446.47,24.5608], Tmin=(100,'K'), Tmax=(740.672,'K')), NASAPolynomial(coeffs=[7.34978,0.0448264,-2.5038e-05,5.10997e-09,-3.6724e-13,-7182.71,-6.02181], Tmin=(740.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.6891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFO) + group(Cds-CdsCsOs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C]1[CH]C2(F)OOC12(943)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {5,S} {8,S} {10,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (24.0771,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,3150,900,1100,212,367,445,1450,792.942,792.942,792.942,792.942,792.943,792.943,792.943,792.943,792.943,792.943],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01799,0.0472587,-3.27286e-05,9.66921e-09,-1.11055e-12,2962.68,18.0674], Tmin=(100,'K'), Tmax=(1964.1,'K')), NASAPolynomial(coeffs=[16.3603,0.0180501,-1.04219e-05,2.09781e-09,-1.46839e-13,-2671.3,-60.8074], Tmin=(1964.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.0771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FC#CC1[C](F)C2OOC12F(1163)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {13,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {11,T}
11 C u0 p0 c0 {3,S} {10,T}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-38.9034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,2175,525,239,401,1367,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943296,0.0747804,-6.9903e-05,3.36382e-08,-6.85072e-12,-4575.3,23.7978], Tmin=(100,'K'), Tmax=(1127.72,'K')), NASAPolynomial(coeffs=[11.5438,0.0371804,-1.98902e-05,4.07217e-09,-2.96294e-13,-6966.15,-28.6173], Tmin=(1127.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.9034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Ct-CtCs) + group(CtCF) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s)"""),
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
    label = 'F[C]=CC1=C(F)C2OO[C]12(1164)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
6  C u1 p0 c0 {4,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {9,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u0 p0 c0 {7,S} {10,D} {12,S}
10 C u1 p0 c0 {2,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (296.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,323,467,575,827,1418,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576476,0.0724549,-6.79416e-05,3.02812e-08,-5.34505e-12,35821.5,25.564], Tmin=(100,'K'), Tmax=(1351.93,'K')), NASAPolynomial(coeffs=[17.8821,0.0212511,-1.11283e-05,2.26461e-09,-1.64078e-13,31142.4,-63.1432], Tmin=(1351.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOO) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[C]=CC1[C](F)C2=C1OO2(1165)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {8,S}
8  C u1 p0 c0 {1,S} {5,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u1 p0 c0 {2,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (417.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,346,659,817,1284,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16115,0.0616886,-4.8894e-05,1.80727e-08,-2.67454e-12,50272.7,25.2798], Tmin=(100,'K'), Tmax=(1559.75,'K')), NASAPolynomial(coeffs=[15.9399,0.0237882,-1.24455e-05,2.49393e-09,-1.77536e-13,45662.4,-52.5887], Tmin=(1559.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ene_m) + radical(CsCdCsF1s) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]C=C1C(F)C2OOC12F(1166)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {9,D} {11,S} {14,S}
11 C u2 p0 c0 {3,S} {10,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-106.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04053,0.0761057,-6.35991e-05,2.56319e-08,-4.36217e-12,-12758.4,20.6606], Tmin=(100,'K'), Tmax=(1301.13,'K')), NASAPolynomial(coeffs=[12.5571,0.0407009,-2.27828e-05,4.71862e-09,-3.43878e-13,-15755.3,-37.9319], Tmin=(1301.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ane) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=CC1C(F)[C]2OOC21F(1167)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {9,S} {13,S}
9  C u1 p0 c0 {5,S} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-49.6116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,316,385,515,654,689,1295,259,529,569,1128,1321,1390,3140,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762666,0.0765436,-6.35673e-05,2.54085e-08,-4.18318e-12,-5855.53,23.7363], Tmin=(100,'K'), Tmax=(1380.01,'K')), NASAPolynomial(coeffs=[14.9277,0.0354857,-1.89393e-05,3.84916e-09,-2.77519e-13,-9765.1,-49.1642], Tmin=(1380.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.6116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=[C]C1[C](F)C2OOC12F(1168)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {13,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {11,D} {14,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-58.0371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74328,0.0780829,-6.83995e-05,2.95073e-08,-5.2896e-12,-6868.8,25.7004], Tmin=(100,'K'), Tmax=(1274.92,'K')), NASAPolynomial(coeffs=[13.9486,0.036652,-1.96542e-05,4.01812e-09,-2.91411e-13,-10236,-41.215], Tmin=(1274.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.0371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Cds_S)"""),
)

species(
    label = 'F[C]=[C]C1C(F)C2OOC12F(1169)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {8,S} {14,S}
10 C u1 p0 c0 {7,S} {11,D}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (1.60673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,250,417,511,1155,1315,1456,3119,1685,370,167,640,1190,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344241,0.0808441,-7.08307e-05,2.96243e-08,-4.98837e-12,323.94,24.9192], Tmin=(100,'K'), Tmax=(1387.74,'K')), NASAPolynomial(coeffs=[17.9811,0.0300076,-1.58817e-05,3.22684e-09,-2.32883e-13,-4571.12,-65.9477], Tmin=(1387.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.60673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[CH]C=C1[C](F)C2OOC12F(1078)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {8,D} {11,S} {13,S}
11 C u1 p0 c0 {3,S} {10,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-208.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,346,659,817,1284,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.319206,0.0797307,-6.44137e-05,2.37613e-08,-3.48069e-12,-24914.2,23.9825], Tmin=(100,'K'), Tmax=(1582.82,'K')), NASAPolynomial(coeffs=[20.7458,0.02811,-1.5494e-05,3.15691e-09,-2.26312e-13,-31380.6,-83.9441], Tmin=(1582.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ane) + radical(CsCdCsF1s) + radical(Csj(Cd-CdH)(F1s)(H))"""),
)

species(
    label = 'FC=CC1[C](F)[C]2OOC21F(1170)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
8  C u1 p0 c0 {5,S} {7,S} {9,S}
9  C u1 p0 c0 {2,S} {6,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {13,S}
11 C u0 p0 c0 {3,S} {10,D} {14,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-109.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07774,0.0745853,-6.32243e-05,2.72394e-08,-5.07876e-12,-13044,24.8334], Tmin=(100,'K'), Tmax=(1187.92,'K')), NASAPolynomial(coeffs=[10.5047,0.0428429,-2.31432e-05,4.7459e-09,-3.45014e-13,-15283.7,-22.2696], Tmin=(1187.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = 'F[C]=CC1[C]2OOC2C1(F)F(1171)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u1 p0 c0 {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-37.7243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,222,329,445,522,589,1214,1475,3010,987.5,1337.5,450,1655,167,640,1190,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823877,0.0762708,-6.60575e-05,2.85808e-08,-5.17536e-12,-4428.54,25.5789], Tmin=(100,'K'), Tmax=(1257.97,'K')), NASAPolynomial(coeffs=[13.062,0.0373569,-1.96567e-05,3.99056e-09,-2.88467e-13,-7507.58,-36.2716], Tmin=(1257.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.7243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]1C(C=C(F)F)[C]2OOC12(1172)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {13,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u1 p0 c0 {1,S} {6,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-66.4388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26633,0.0694895,-5.45693e-05,2.14684e-08,-3.64337e-12,-7900.66,26.7559], Tmin=(100,'K'), Tmax=(1283.77,'K')), NASAPolynomial(coeffs=[10.3055,0.0413248,-2.16606e-05,4.37876e-09,-3.15338e-13,-10221.5,-19.111], Tmin=(1283.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.4388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFF) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
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
    E0 = (-29.4649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (276.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (130.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (130.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (130.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-21.1806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (33.9352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-4.49169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (37.6947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (71.9447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-29.4649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (127.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (110.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (173.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (141.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (189.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (390.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (264.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (334.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (104.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (90.5891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (140.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (150.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (155.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (14.8329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (202.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (201.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['C2HF(58)', 'FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]=CC1[C](F)C(F)C12OO2(1153)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['F[C]C=CC1(F)C2OOC21F(1154)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['F[C]=CC1(F)[CH]C2(F)OOC12(1019)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['F[C]=CC1C2(F)[CH]OOC12F(1155)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['FC1=CC2C3(F)OOC3C12F(1156)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_noH;Cdsinglepri_rad_out]
Euclidian distance = 1.7320508075688772
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['FC=CC1=C(F)C2OOC12F(1157)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['FC=CC1C(F)=C2OOC21F(1021)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[C]=CC1[C](F)OOC=C1F(1158)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c6_gamma;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C]=CC=C(F)C1OO[C]1F(1159)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC1(F)C=C(F)C1C=[C]F(1160)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.10864e+11,'s^-1'), n=0.310877, Ea=(90.4303,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 88.2 to 90.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]F(252)', 'FC1=CC2(F)OOC12(939)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(9.66596,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'F[C]=CC1=C(F)C2OOC12F(1161)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(205.9,'m^3/(mol*s)'), n=1.596, Ea=(16.9876,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'F[C]=CC1C(F)=C2OOC21F(1162)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.296163,'m^3/(mol*s)'), n=2.44476, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_Ext-9R!H-R_Ext-8R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_Ext-9R!H-R_Ext-8R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C2HF(58)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.0408,'m^3/(mol*s)'), n=2.24, Ea=(5.77665,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'FC#CC1[C](F)C2OOC12F(1163)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.29124e+14,'m^3/(mol*s)'), n=-1.78703, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]F(252)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'F[C]=CC1=C(F)C2OO[C]12(1164)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(232.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'F[C]=CC1[C](F)C2=C1OO2(1165)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(181.78,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['F[C]C=C1C(F)C2OOC12F(1166)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]=CC1C(F)[C]2OOC21F(1167)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.26888e+09,'s^-1'), n=0.972714, Ea=(123.905,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_NDMustO;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['FC=[C]C1[C](F)C2OOC12F(1168)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]=[C]C1C(F)C2OOC12F(1169)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.1e+10,'s^-1'), n=0.78, Ea=(132.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_23cy4;Cd_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['F[CH]C=C1[C](F)C2OOC12F(1078)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(323312,'s^-1'), n=2.11016, Ea=(185.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['FC=CC1[C](F)[C]2OOC21F(1170)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_single;XH_out]
Euclidian distance = 2.8284271247461903
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]=CC1[C]2OOC2C1(F)F(1171)'],
    products = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(223.716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[C]=CC1[C](F)C2OOC12F(1020)'],
    products = ['F[C]1C(C=C(F)F)[C]2OOC12(1172)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(231.347,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #222',
    isomers = [
        'F[C]=CC1[C](F)C2OOC12F(1020)',
    ],
    reactants = [
        ('C2HF(58)', 'FC1=CC2(F)OOC12(939)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #222',
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

