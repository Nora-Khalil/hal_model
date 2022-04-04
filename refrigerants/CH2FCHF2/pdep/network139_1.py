species(
    label = 'F[C]=CC1[C](F)C=C1F(601)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-4.43144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,180,1017.95,1021.55,1024.56,1028.97],'cm^-1')),
        HinderedRotor(inertia=(0.00875121,'amu*angstrom^2'), symmetry=1, barrier=(6.49343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03883,0.0653693,-6.3397e-05,3.16221e-08,-6.40058e-12,-426.579,24.2569], Tmin=(100,'K'), Tmax=(1176.62,'K')), NASAPolynomial(coeffs=[13.0006,0.0247048,-1.15565e-05,2.24966e-09,-1.59774e-13,-3241.48,-35.3972], Tmin=(1176.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.43144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'F[C]C=CC1(F)C=C1F(810)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {5,D}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u2 p0 c0 {3,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (73.3528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,323,467,575,827,1418,2995,3025,975,1000,1300,1375,400,500,1630,1680,301.889,303.362,307.159,307.601,321.954],'cm^-1')),
        HinderedRotor(inertia=(0.00606659,'amu*angstrom^2'), symmetry=1, barrier=(8.90411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335337,'amu*angstrom^2'), symmetry=1, barrier=(22.2173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8862,0.0749346,-0.000100057,7.46876e-08,-2.30323e-11,8928.61,25.903], Tmin=(100,'K'), Tmax=(783.555,'K')), NASAPolynomial(coeffs=[9.48359,0.0310485,-1.60495e-05,3.21703e-09,-2.30623e-13,7581.21,-13.478], Tmin=(783.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.3528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs(F)-Cd) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C]=CC1(F)[CH]C(F)=C1(600)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-8.8628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26147,0.0624005,-5.85646e-05,2.86969e-08,-5.8112e-12,-969.192,22.4732], Tmin=(100,'K'), Tmax=(1160.37,'K')), NASAPolynomial(coeffs=[11.4637,0.027232,-1.31029e-05,2.57796e-09,-1.83943e-13,-3336.87,-28.2639], Tmin=(1160.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.8628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC1=CC2(F)C(F)=CC12(602)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-241.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05505,0.05333,-2.54873e-05,-7.57033e-09,6.73605e-12,-28886,22.1434], Tmin=(100,'K'), Tmax=(1065.09,'K')), NASAPolynomial(coeffs=[16.3701,0.0184805,-8.33028e-06,1.67206e-09,-1.23442e-13,-33434.1,-58.7441], Tmin=(1065.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + polycyclic(s2_4_4_diene_1_4)"""),
)

species(
    label = 'FC=CC1=C(F)C=C1F(811)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,D} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {3,S} {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-101.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30669,0.0582106,-4.819e-05,2.00995e-08,-3.42573e-12,-12064.3,22.2422], Tmin=(100,'K'), Tmax=(1366.43,'K')), NASAPolynomial(coeffs=[12.7011,0.0248552,-1.1574e-05,2.23494e-09,-1.57248e-13,-15178.2,-36.2868], Tmin=(1366.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(CdCCF) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cd-Cd-Cd(F))"""),
)

species(
    label = 'FC=CC1C(F)=C=C1F(604)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {4,S} {8,D} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {9,D}
7  C u0 p0 c0 {2,S} {4,S} {9,D}
8  C u0 p0 c0 {3,S} {5,D} {12,S}
9  C u0 p0 c0 {6,D} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (9.65833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.893644,0.0766023,-0.000120672,1.10388e-07,-4.02334e-11,1265.46,23.9617], Tmin=(100,'K'), Tmax=(805.618,'K')), NASAPolynomial(coeffs=[6.2286,0.0358287,-1.81571e-05,3.54497e-09,-2.47656e-13,869.434,2.25377], Tmin=(805.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.65833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(CdCddCF) + group(CdCddCF) + group(CdCFH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'F[C]=CC1C2(F)[CH]C12F(812)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {11,S}
8  C u0 p0 c0 {6,S} {9,D} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (124.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2745,0.0658883,-7.43691e-05,4.7225e-08,-1.26873e-11,15117.7,25.5382], Tmin=(100,'K'), Tmax=(883.005,'K')), NASAPolynomial(coeffs=[8.84525,0.0315927,-1.61092e-05,3.23854e-09,-2.3357e-13,13780.7,-10.0443], Tmin=(883.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCCF) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_3_3_ane) + radical(bicyclo[1.1.0]butane-secondary) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1C2C=C(F)C1[C]2F(813)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {11,S}
6  C u1 p0 c0 {1,S} {4,S} {5,S}
7  C u1 p0 c0 {2,S} {4,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-13.8556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27729,0.0530371,-3.72864e-05,1.21375e-08,-1.55331e-12,-1562.82,22.8521], Tmin=(100,'K'), Tmax=(1824.83,'K')), NASAPolynomial(coeffs=[17.299,0.0179178,-8.41854e-06,1.59114e-09,-1.08466e-13,-7410.18,-64.0799], Tmin=(1824.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.8556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s3_4_5_ene_1) + radical(CsCsCsF1s) + radical(CsCsCsF1s)"""),
)

species(
    label = 'F[C]1[CH]C2(F)C(F)=CC12(814)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
6  C u1 p0 c0 {2,S} {4,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-12.5277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08149,0.0515182,-1.94902e-05,-6.13856e-08,7.24498e-11,-1447.27,21.0143], Tmin=(100,'K'), Tmax=(444.452,'K')), NASAPolynomial(coeffs=[4.83874,0.0392947,-2.07317e-05,4.21825e-09,-3.05841e-13,-1816.72,8.54904], Tmin=(444.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.5277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(CsCsCsF1s) + radical(cyclobutane) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C]=CC=C(F)C=[C]F(815)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,D} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {4,D} {6,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u1 p0 c0 {2,S} {6,D}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (78.4413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,280,518,736,852,873,125,209,569,711,1121,1259,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.737079,'amu*angstrom^2'), symmetry=1, barrier=(16.9469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732099,'amu*angstrom^2'), symmetry=1, barrier=(16.8324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.150978,0.0896388,-0.000129987,9.64499e-08,-2.82206e-11,9568.42,26.5934], Tmin=(100,'K'), Tmax=(839.746,'K')), NASAPolynomial(coeffs=[14.1063,0.023162,-1.12372e-05,2.1715e-09,-1.51837e-13,7224.74,-38.2948], Tmin=(839.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.4413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCFH) + radical(Cdj(Cd-CdH)(F1s)) + radical(Cdj(Cd-CdH)(F1s))"""),
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
    label = 'F[C]=CC1=C(F)C=C1F(816)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,D} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (154.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([217,343,449,587,660,812,790,914,798,948,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,1441.05,1441.06,1441.27,1441.79],'cm^-1')),
        HinderedRotor(inertia=(0.383554,'amu*angstrom^2'), symmetry=1, barrier=(8.81867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31281,0.0644357,-7.69352e-05,5.17056e-08,-1.45411e-11,18685.2,22.1684], Tmin=(100,'K'), Tmax=(851.017,'K')), NASAPolynomial(coeffs=[8.82894,0.0291068,-1.46623e-05,2.92087e-09,-2.09305e-13,17406,-12.8798], Tmin=(851.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(CdCCF) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cd-Cd-Cd(F)) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[C]=CC1C(F)=C=C1F(817)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {8,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u0 p0 c0 {5,D} {6,D}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (259.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,104,186,245,407,330,466,761,907,1218,1388,3010,987.5,1337.5,450,1655,167,640,1190,180,180,1906.72,1906.85,1906.93,1907.27],'cm^-1')),
        HinderedRotor(inertia=(0.320365,'amu*angstrom^2'), symmetry=1, barrier=(7.36582,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844348,0.079822,-0.000138947,1.31146e-07,-4.76844e-11,31347.4,24.6291], Tmin=(100,'K'), Tmax=(834.906,'K')), NASAPolynomial(coeffs=[6.02844,0.0335795,-1.74093e-05,3.39028e-09,-2.34889e-13,31227.9,5.02221], Tmin=(834.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(CdCddCF) + group(CdCddCF) + group(CdCFH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'FC#CC1[C](F)C=C1F(818)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {3,S} {8,T}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-20.9734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,2175,525,239,401,1367,180,180,757.168,757.171,757.2,757.21,757.224,757.401],'cm^-1')),
        HinderedRotor(inertia=(0.0114813,'amu*angstrom^2'), symmetry=1, barrier=(4.6728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1821,0.063178,-6.38522e-05,3.27305e-08,-6.78648e-12,-2422.05,22.054], Tmin=(100,'K'), Tmax=(1150.88,'K')), NASAPolynomial(coeffs=[12.8076,0.0227721,-1.11891e-05,2.22438e-09,-1.59776e-13,-5097.97,-35.6662], Tmin=(1150.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.9734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(CtCF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s)"""),
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
    label = 'F[C]=CC1=C(F)C=[C]1(819)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,D} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,D} {5,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u0 p0 c0 {3,S} {8,D} {10,S}
7  C u1 p0 c0 {3,S} {5,D}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (569.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,180,180,1529.67,1533.71,1535.31],'cm^-1')),
        HinderedRotor(inertia=(0.00186076,'amu*angstrom^2'), symmetry=1, barrier=(3.09908,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51224,0.0585642,-7.5261e-05,5.42694e-08,-1.60413e-11,68534.8,22.0406], Tmin=(100,'K'), Tmax=(819.94,'K')), NASAPolynomial(coeffs=[8.69436,0.0235268,-1.11634e-05,2.15359e-09,-1.51095e-13,67357,-11.1832], Tmin=(819.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cd-Cd-Cd(F)) + radical(cyclobutadiene-C1) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[C]=CC1[C]=C=C1F(820)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {7,D}
5  C u0 p0 c0 {3,S} {8,D} {10,S}
6  C u1 p0 c0 {3,S} {7,D}
7  C u0 p0 c0 {4,D} {6,D}
8  C u1 p0 c0 {2,S} {5,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (688.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,145,326,398,834,1303,3010,987.5,1337.5,450,1655,167,640,1190,180.022,180.025,180.201,180.216,198.558,1927.67,1927.99,1930.69],'cm^-1')),
        HinderedRotor(inertia=(0.000677461,'amu*angstrom^2'), symmetry=1, barrier=(1.79032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2877,0.0704,-0.000127371,1.23805e-07,-4.56147e-11,82903.9,23.8502], Tmin=(100,'K'), Tmax=(847.917,'K')), NASAPolynomial(coeffs=[4.43844,0.0315772,-1.63064e-05,3.15654e-09,-2.17347e-13,83230.9,14.2484], Tmin=(847.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsH) + group(CdCddCF) + group(Cds-CdsCsH) + group(CdCFH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]C=C1C(F)=CC1F(821)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {8,D}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {5,D} {9,S} {12,S}
9  C u2 p0 c0 {3,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-18.8537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26339,0.0635314,-5.82187e-05,2.65895e-08,-4.99245e-12,-2172.03,20.2743], Tmin=(100,'K'), Tmax=(1239.6,'K')), NASAPolynomial(coeffs=[12.5609,0.0270754,-1.41036e-05,2.86353e-09,-2.07353e-13,-4972.86,-36.656], Tmin=(1239.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCCF) + ring(Cyclobutene) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=CC1C(F)=[C]C1F(822)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u1 p0 c0 {5,S} {6,D}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (122.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,164,312,561,654,898,1207,1299,3167,246,474,533,1155,3010,987.5,1337.5,450,1655,167,640,1190,196.161,196.162,196.186,1002.73,1581.27,1581.28,1581.28],'cm^-1')),
        HinderedRotor(inertia=(0.00438124,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1602,0.0662397,-7.03181e-05,4.0141e-08,-9.47372e-12,14822.6,23.5591], Tmin=(100,'K'), Tmax=(1008.61,'K')), NASAPolynomial(coeffs=[10.7641,0.0281515,-1.36727e-05,2.69922e-09,-1.9304e-13,12885.4,-22.8563], Tmin=(1008.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=[C]C1[C](F)C=C1F(823)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {3,S} {9,D} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-16.7076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11692,0.0659591,-6.6875e-05,3.58139e-08,-7.88237e-12,-1907.73,23.9644], Tmin=(100,'K'), Tmax=(1079.93,'K')), NASAPolynomial(coeffs=[11.5911,0.0271636,-1.29893e-05,2.54928e-09,-1.81806e-13,-4170.03,-27.3731], Tmin=(1079.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.7076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cds_S)"""),
)

species(
    label = 'F[C]=[C]C1C(F)=CC1F(824)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {12,S}
8  C u1 p0 c0 {4,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (95.4634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,284,328,853,1146,1135,1297,3239,323,467,575,827,1418,1685,370,167,640,1190,283.797,283.8,283.814,283.814,1596,1596,1596,1596.01],'cm^-1')),
        HinderedRotor(inertia=(0.146032,'amu*angstrom^2'), symmetry=1, barrier=(8.34622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.016,0.0711932,-8.90105e-05,6.28465e-08,-1.8418e-11,11584.1,23.9536], Tmin=(100,'K'), Tmax=(821.72,'K')), NASAPolynomial(coeffs=[9.31592,0.0307911,-1.526e-05,3.01313e-09,-2.14571e-13,10220.1,-14.459], Tmin=(821.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.4634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=C[C]1[C](F)C=C1F(762)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u1 p0 c0 {5,S} {6,S} {8,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {3,S} {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-122.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([346,659,817,1284,271,519,563,612,1379,2950,1000,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,786.694,786.702,786.751,786.765,786.789,786.815],'cm^-1')),
        HinderedRotor(inertia=(0.133061,'amu*angstrom^2'), symmetry=1, barrier=(3.05934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10209,0.0598507,-5.03206e-05,2.1468e-08,-3.70273e-12,-14634.2,21.7916], Tmin=(100,'K'), Tmax=(1370.29,'K')), NASAPolynomial(coeffs=[13.6943,0.0230927,-1.0083e-05,1.89174e-09,-1.31151e-13,-18085.2,-42.9252], Tmin=(1370.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Allyl_T) + radical(CsCdCsF1s)"""),
)

species(
    label = 'FC=CC1[C](F)[C]=C1F(825)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {4,S} {8,D} {11,S}
6  C u1 p0 c0 {1,S} {4,S} {9,S}
7  C u0 p0 c0 {2,S} {4,S} {9,D}
8  C u0 p0 c0 {3,S} {5,D} {12,S}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (10.2471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3010,987.5,1337.5,450,1655,346,659,817,1284,246,474,533,1155,194,682,905,1196,1383,3221,180,180,180,907.519,909.202,909.22,911.887,913.685],'cm^-1')),
        HinderedRotor(inertia=(0.104819,'amu*angstrom^2'), symmetry=1, barrier=(2.40999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982255,0.064303,-5.97028e-05,2.8004e-08,-5.28428e-12,1342.72,23.8755], Tmin=(100,'K'), Tmax=(1263.61,'K')), NASAPolynomial(coeffs=[14.0935,0.0228001,-1.0437e-05,2.01258e-09,-1.42123e-13,-1970.86,-42.4466], Tmin=(1263.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.2471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
)

species(
    label = 'F[C]=CC1[C]=CC1(F)F(826)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u1 p0 c0 {4,S} {6,D}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (90.9739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,274,345,380,539,705,1166,1213,3010,987.5,1337.5,450,1655,167,640,1190,401.358,401.361,401.365,401.37,401.37,401.377,1622.37,1622.38,1622.38,1622.38],'cm^-1')),
        HinderedRotor(inertia=(0.0495997,'amu*angstrom^2'), symmetry=1, barrier=(5.66999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13126,0.0673051,-7.39864e-05,4.41779e-08,-1.09125e-11,11041.4,23.6755], Tmin=(100,'K'), Tmax=(966.564,'K')), NASAPolynomial(coeffs=[10.4842,0.0285987,-1.39178e-05,2.74652e-09,-1.96279e-13,9233.4,-21.129], Tmin=(966.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.9739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(cyclobutene-vinyl) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]1C=[C]C1C=C(F)F(827)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (5.40367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,339.744,340.017,340.676,341.301,341.527,876.474,1435.63,1435.93,1436.06,2223.72],'cm^-1')),
        HinderedRotor(inertia=(0.539392,'amu*angstrom^2'), symmetry=1, barrier=(44.8135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23562,0.0615115,-5.56891e-05,2.60395e-08,-5.00121e-12,748.744,24.4413], Tmin=(100,'K'), Tmax=(1223.5,'K')), NASAPolynomial(coeffs=[12.0623,0.0261159,-1.22946e-05,2.39465e-09,-1.69837e-13,-1900.57,-29.9754], Tmin=(1223.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.40367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(cyclobutene-vinyl)"""),
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
    E0 = (37.6719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (139.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (71.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-82.8478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-12.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-66.1588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (122.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-11.8871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-89.8769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (129.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (292.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (279.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (384.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (98.9116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (104.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (350.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (343.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (512.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (42.3676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (218.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (78.6149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (141.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (94.1578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (102.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (186.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (156.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['C2HF(58)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(128.804,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 128.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]C=CC1(F)C=C1F(810)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.65366e+10,'s^-1'), n=0.611, Ea=(152.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCsCJ;CsJ-CdH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CdH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.13382e+11,'s^-1'), n=0.663, Ea=(162.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-OneDe;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['FC1=CC2(F)C(F)=CC12(602)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_noH;Cdsinglepri_rad_out]
Euclidian distance = 1.7320508075688772
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['FC=CC1=C(F)C=C1F(811)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['FC=CC1C(F)=C=C1F(604)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['F[C]=CC1C2(F)[CH]C12F(812)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['F[C]1C2C=C(F)C1[C]2F(813)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.18, Ea=(79.245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra;radadd_intra_cdsingle] for rate rule [Rn2c4_beta;doublebond_intra;radadd_intra_cdsingle]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['F[C]1[CH]C2(F)C(F)=CC12(814)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.2e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_cdsingle]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C]=CC=C(F)C=[C]F(815)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingle]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]F(252)', 'FC1=CC(F)=C1(307)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'F[C]=CC1=C(F)C=C1F(816)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'F[C]=CC1C(F)=C=C1F(817)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.10362,'m^3/(mol*s)'), n=2.13572, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R_Ext-5R!H-R_Ext-2CS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R_Ext-5R!H-R_Ext-2CS-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C2HF(58)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0816,'m^3/(mol*s)'), n=2.24, Ea=(3.47577,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'FC#CC1[C](F)C=C1F(818)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.29124e+14,'m^3/(mol*s)'), n=-1.78703, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]F(252)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'F[C]=CC1=C(F)C=[C]1(819)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(142.029,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'F[C]=CC1[C]=C=C1F(820)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(192.236,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['F[C]C=C1C(F)=CC1F(821)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.41827e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C]=CC1C(F)=[C]C1F(822)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['FC=[C]C1[C](F)C=C1F(823)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]=[C]C1C(F)=CC1F(824)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.1e+10,'s^-1'), n=0.78, Ea=(132.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_23cy4;Cd_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]=CC1[C](F)C=C1F(601)'],
    products = ['FC=C[C]1[C](F)C=C1F(762)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(323312,'s^-1'), n=2.11016, Ea=(185.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FC=CC1[C](F)[C]=C1F(825)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.02951e+08,'s^-1'), n=1.2385, Ea=(178.845,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_single] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cd_H_out_single]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C]=CC1[C]=CC1(F)F(826)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(182.661,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]1C=[C]C1C=C(F)F(827)'],
    products = ['F[C]=CC1[C](F)C=C1F(601)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(237.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #139',
    isomers = [
        'F[C]=CC1[C](F)C=C1F(601)',
    ],
    reactants = [
        ('C2HF(58)', 'FC1=CC(F)=C1(307)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #139',
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

