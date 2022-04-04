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
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,3010,987.5,1337.5,450,1655,167,640,1190,564.012,564.333,564.399,564.745,564.813,564.9],'cm^-1')),
        HinderedRotor(inertia=(0.000529874,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26147,0.0624005,-5.85646e-05,2.86969e-08,-5.8112e-12,-969.192,22.4732], Tmin=(100,'K'), Tmax=(1160.37,'K')), NASAPolynomial(coeffs=[11.4637,0.027232,-1.31029e-05,2.57796e-09,-1.83943e-13,-3336.87,-28.2639], Tmin=(1160.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.8628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'F[C]C=C(F)C1C=C1F(828)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u2 p0 c0 {3,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (132.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,260,386,409,525,515,635,761,893,1354,1482,3010,987.5,1337.5,450,1655,180,180,180,437.492,730.221,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.17901,'amu*angstrom^2'), symmetry=1, barrier=(4.1158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17901,'amu*angstrom^2'), symmetry=1, barrier=(4.1158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753059,0.0780898,-0.000117201,1.00709e-07,-3.52918e-11,16073.1,25.2563], Tmin=(100,'K'), Tmax=(760.093,'K')), NASAPolynomial(coeffs=[8.43759,0.0323629,-1.6528e-05,3.25848e-09,-2.29841e-13,15057.7,-8.70446], Tmin=(760.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CsCFHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cd(F)-Cd) + radical(CsCCl_triplet)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3493.2,'J/mol'), sigma=(5.43195,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.63 K, Pc=49.45 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03883,0.0653693,-6.3397e-05,3.16221e-08,-6.40058e-12,-426.579,24.2569], Tmin=(100,'K'), Tmax=(1176.62,'K')), NASAPolynomial(coeffs=[13.0006,0.0247048,-1.15565e-05,2.24966e-09,-1.59774e-13,-3241.48,-35.3972], Tmin=(1176.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.43144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=CC1=CC(F)(F)[CH]1(829)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
5  C u0 p0 c0 {6,S} {7,D} {8,S}
6  C u1 p0 c0 {4,S} {5,S} {11,S}
7  C u0 p0 c0 {4,S} {5,D} {10,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-12.8712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2750,3150,900,1100,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,840.437,840.442,840.477,840.491,840.498,840.5],'cm^-1')),
        HinderedRotor(inertia=(0.0100721,'amu*angstrom^2'), symmetry=1, barrier=(5.04858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597956,0.0655763,-5.6234e-05,2.06177e-08,-2.04981e-12,-1417.67,22.5342], Tmin=(100,'K'), Tmax=(1093.42,'K')), NASAPolynomial(coeffs=[17.4365,0.0168203,-6.96752e-06,1.32204e-09,-9.42375e-14,-5867.73,-63.7164], Tmin=(1093.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.8712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFF) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(cyclobutene-allyl) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'FC1=CC2(F)C=C(F)C12(603)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-256.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919095,0.0576298,-3.89904e-05,7.33963e-09,1.28198e-12,-30702.3,20.5662], Tmin=(100,'K'), Tmax=(1143.28,'K')), NASAPolynomial(coeffs=[16.5101,0.0188224,-8.72683e-06,1.73518e-09,-1.25898e-13,-35296,-61.2381], Tmin=(1143.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_diene_1_4)"""),
)

species(
    label = 'F[C]=CC1(F)C2[C](F)C21(830)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u1 p0 c0 {2,S} {4,S} {5,S}
8  C u0 p0 c0 {6,S} {9,D} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (73.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08122,0.0709785,-9.00848e-05,6.57932e-08,-2.02754e-11,8971.26,23.5057], Tmin=(100,'K'), Tmax=(777.393,'K')), NASAPolynomial(coeffs=[8.3935,0.0333545,-1.74893e-05,3.53867e-09,-2.55421e-13,7834.34,-9.93067], Tmin=(777.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + polycyclic(s2_3_3_ane) + radical(CsCsCsF1s) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC1=CC2(F)[CH]C1(F)[CH]2(831)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {4,S} {5,S} {11,S}
7  C u1 p0 c0 {4,S} {5,S} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (25.5373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52404,0.0606851,-7.63557e-05,6.88411e-08,-2.82646e-11,3154.74,20.7259], Tmin=(100,'K'), Tmax=(656.183,'K')), NASAPolynomial(coeffs=[4.17537,0.0406691,-2.17901e-05,4.45304e-09,-3.2317e-13,2889.76,9.68401], Tmin=(656.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.5373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s3_4_5_ene_1) + radical(bicyclo[2.1.1]hex-2-ene-C5) + radical(bicyclo[2.1.1]hex-2-ene-C5) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C]1[CH]C2(F)C=C(F)C12(832)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {7,S} {11,S}
7  C u1 p0 c0 {2,S} {5,S} {6,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-27.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26748,0.067008,-9.41488e-05,8.72608e-08,-3.45568e-11,-3234.84,22.3946], Tmin=(100,'K'), Tmax=(714.522,'K')), NASAPolynomial(coeffs=[4.7543,0.0398745,-2.12039e-05,4.28834e-09,-3.08181e-13,-3538.76,8.10482], Tmin=(714.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(cyclobutane) + radical(CsCsCsF1s)"""),
)

species(
    label = '[CH]=C(F)C=C(F)C=[C]F(833)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,D} {7,S}
5  C u0 p0 c0 {4,D} {6,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u1 p0 c0 {6,D} {12,S}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (82.5719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2995,3025,975,1000,1300,1375,400,500,1630,1680,250,446,589,854,899,3120,650,792.5,1650,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.785708,'amu*angstrom^2'), symmetry=1, barrier=(18.065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785635,'amu*angstrom^2'), symmetry=1, barrier=(18.0633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291007,0.0865185,-0.000117049,8.07282e-08,-2.20993e-11,10060.4,26.186], Tmin=(100,'K'), Tmax=(892.902,'K')), NASAPolynomial(coeffs=[14.2258,0.0240909,-1.21712e-05,2.41939e-09,-1.72903e-13,7572.05,-39.462], Tmin=(892.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.5719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCFH) + radical(Cdj(Cd-CdF1s)(H)) + radical(Cdj(Cd-CdH)(F1s))"""),
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
    label = 'F[C]=CC1=CC(F)=C1(834)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,D} {5,S} {7,S}
4  C u0 p0 c0 {3,D} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (349.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,280,518,736,852,873,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,180,368.755,1417.62,1419.12,1419.75,1420.13],'cm^-1')),
        HinderedRotor(inertia=(0.00223864,'amu*angstrom^2'), symmetry=1, barrier=(3.19402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5703,0.0537691,-5.01706e-05,2.49382e-08,-5.10057e-12,42144.3,20.6579], Tmin=(100,'K'), Tmax=(1157.72,'K')), NASAPolynomial(coeffs=[10.4881,0.0229583,-1.02512e-05,1.95135e-09,-1.36841e-13,40079.4,-23.6714], Tmin=(1157.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cd-Cd-Cd(F)) + radical(Cdj(Cd-CdH)(F1s))"""),
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
    label = 'F[C]=CC1(F)C=C=C1(835)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {10,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {3,S} {8,D} {9,S}
7  C u0 p0 c0 {4,D} {5,D}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (403.626,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,180,180,1503.49,1503.97,1503.97],'cm^-1')),
        HinderedRotor(inertia=(0.124376,'amu*angstrom^2'), symmetry=1, barrier=(2.85966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12672,0.0738503,-0.000129026,1.26352e-07,-4.77531e-11,48638.1,21.7045], Tmin=(100,'K'), Tmax=(820.952,'K')), NASAPolynomial(coeffs=[4.22852,0.0362408,-1.92035e-05,3.79028e-09,-2.65229e-13,48886.9,11.9692], Tmin=(820.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'FC#CC1(F)[CH]C(F)=C1(836)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {3,S} {8,T}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1.00528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,2175,525,239,401,1367,611.668,613.017,614.378,614.699,614.841,615.228],'cm^-1')),
        HinderedRotor(inertia=(0.000448746,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.075,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39321,0.0606485,-6.41417e-05,3.69744e-08,-8.84425e-12,-29.7768,20.6009], Tmin=(100,'K'), Tmax=(995.044,'K')), NASAPolynomial(coeffs=[9.89109,0.0264877,-1.26454e-05,2.47248e-09,-1.75814e-13,-1720.93,-20.3542], Tmin=(995.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.00528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Ct-CtCs) + group(CtCF) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl)"""),
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
    label = 'F[C]=CC1=C=C(F)[CH]1(837)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {6,S} {7,D}
4  C u1 p0 c0 {3,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {3,S} {8,D} {10,S}
7  C u0 p0 c0 {3,D} {5,D}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (532.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,132,421,945,1169,1268,3010,987.5,1337.5,450,1655,167,640,1190,180,180,660.296,660.307,660.526,660.637,660.637,660.66],'cm^-1')),
        HinderedRotor(inertia=(0.0111676,'amu*angstrom^2'), symmetry=1, barrier=(3.45482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30086,0.056041,-5.67469e-05,2.90599e-08,-5.88353e-12,64091.2,21.7437], Tmin=(100,'K'), Tmax=(1199.77,'K')), NASAPolynomial(coeffs=[13.3146,0.0159879,-6.67136e-06,1.23504e-09,-8.5618e-14,61208.4,-38.4037], Tmin=(1199.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCddCF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(C=CCJC=C) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[C]=CC1(F)[C]=C=C1(838)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {3,S} {8,D} {10,S}
6  C u1 p0 c0 {3,S} {7,D}
7  C u0 p0 c0 {4,D} {6,D}
8  C u1 p0 c0 {2,S} {5,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (641.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,180,180,1744.46,1748.65],'cm^-1')),
        HinderedRotor(inertia=(0.156737,'amu*angstrom^2'), symmetry=1, barrier=(3.60369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07441,0.0786708,-0.000154517,1.56342e-07,-5.88853e-11,77242.3,23.0593], Tmin=(100,'K'), Tmax=(847.532,'K')), NASAPolynomial(coeffs=[3.19923,0.035461,-1.93166e-05,3.79986e-09,-2.63251e-13,78073.8,20.1903], Tmin=(847.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=C=C1[CH]C(F)=C1(839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {7,D}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u0 p0 c0 {3,D} {8,D}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (287.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,271,519,563,612,1379,540,610,2055,137,207,812,816.331,816.332,816.332,816.332,816.332,816.333,816.333,816.333,816.333],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54668,0.0369867,1.87236e-05,-6.38981e-08,3.08564e-11,34685,19.4145], Tmin=(100,'K'), Tmax=(927.476,'K')), NASAPolynomial(coeffs=[19.1942,0.00451667,6.58646e-07,-1.81472e-10,7.07037e-15,29534.5,-74.5148], Tmin=(927.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCddFH) + group(Cdd-CdsCds) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Cdj(Cdd-Cd)(F1s))"""),
)

species(
    label = 'FC=[C]C1(F)[CH]C(F)=C1(840)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {3,S} {9,D} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-21.1389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31051,0.0632964,-6.29593e-05,3.38945e-08,-7.65749e-12,-2449.01,22.2876], Tmin=(100,'K'), Tmax=(1040.49,'K')), NASAPolynomial(coeffs=[10.0423,0.0297274,-1.45638e-05,2.88546e-09,-2.06708e-13,-4266.03,-20.1847], Tmin=(1040.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.1389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Cds_S)"""),
)

species(
    label = 'F[C]=CC1(F)[C]=C(F)C1(841)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u1 p0 c0 {4,S} {6,D}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (79.0799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,246,474,533,1155,3010,987.5,1337.5,450,1655,167,640,1190,318.184,318.184,318.184,318.184,818.911,2246.44,2246.44],'cm^-1')),
        HinderedRotor(inertia=(0.453286,'amu*angstrom^2'), symmetry=1, barrier=(32.5653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930077,0.0732572,-9.53521e-05,6.99363e-08,-2.11845e-11,9616.5,25.2539], Tmin=(100,'K'), Tmax=(797.883,'K')), NASAPolynomial(coeffs=[9.42784,0.0306542,-1.52568e-05,3.01073e-09,-2.1403e-13,8260.5,-13.8238], Tmin=(797.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.0799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-vinyl) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=[C]C1(F)C=C(F)C1(842)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u1 p0 c0 {4,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (64.4359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2950,3150,900,1000,1100,323,467,575,827,1418,1685,370,167,640,1190,180,180,180,962.264,1213.34,1216.52,1221.09],'cm^-1')),
        HinderedRotor(inertia=(0.329469,'amu*angstrom^2'), symmetry=1, barrier=(7.57515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930077,0.0732572,-9.53521e-05,6.99363e-08,-2.11845e-11,7855.23,25.2539], Tmin=(100,'K'), Tmax=(797.883,'K')), NASAPolynomial(coeffs=[9.42784,0.0306542,-1.52568e-05,3.01073e-09,-2.1403e-13,6499.23,-13.8238], Tmin=(797.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=CC1(F)[C]=C(F)[CH]1(843)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u1 p0 c0 {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {6,D} {12,S}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-6.49495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,3010,987.5,1337.5,450,1655,271,519,563,612,1379,194,682,905,1196,1383,3221,458.79,459.639,460.223,460.305,461.235],'cm^-1')),
        HinderedRotor(inertia=(0.000798664,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31051,0.0632964,-6.29593e-05,3.38945e-08,-7.65749e-12,-687.741,22.2876], Tmin=(100,'K'), Tmax=(1040.49,'K')), NASAPolynomial(coeffs=[10.0423,0.0297274,-1.45638e-05,2.88546e-09,-2.06708e-13,-2504.76,-20.1847], Tmin=(1040.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.49495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'F[C]C=C1C=C(F)C1F(844)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {8,D}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {5,D} {9,S} {12,S}
9  C u2 p0 c0 {3,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-17.1298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866638,0.0659135,-6.23469e-05,2.8901e-08,-5.32425e-12,-1945.05,22.6027], Tmin=(100,'K'), Tmax=(1300.7,'K')), NASAPolynomial(coeffs=[15.7146,0.0202519,-9.68848e-06,1.91105e-09,-1.36641e-13,-5807.57,-52.9334], Tmin=(1300.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.1298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C](F)C=C1[CH]C(F)=C1(639)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,D}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,D}
8  C u0 p0 c0 {4,D} {9,S} {12,S}
9  C u1 p0 c0 {2,S} {3,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-158.312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,271,519,563,612,1379,3010,987.5,1337.5,450,1655,161,297,490,584,780,1358,689.276,689.278,689.338,689.421,689.619,689.686,689.697,689.76,689.772],'cm^-1')),
        HinderedRotor(inertia=(0.000355053,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.992457,0.0480013,8.87809e-06,-5.76413e-08,2.8524e-11,-18915.6,23.0978], Tmin=(100,'K'), Tmax=(948.008,'K')), NASAPolynomial(coeffs=[20.732,0.00986862,-2.23425e-06,4.1772e-10,-3.68463e-14,-24687.4,-81.7819], Tmin=(948.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]=CC1(F)C=[C]C1F(845)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u1 p0 c0 {5,S} {6,D}
9  C u1 p0 c0 {3,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (113.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,883.076,883.444,883.612],'cm^-1')),
        HinderedRotor(inertia=(0.18102,'amu*angstrom^2'), symmetry=1, barrier=(4.162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17663,0.0675772,-7.66199e-05,4.82141e-08,-1.27031e-11,13742.3,24.6234], Tmin=(100,'K'), Tmax=(904.038,'K')), NASAPolynomial(coeffs=[9.5134,0.0306898,-1.54144e-05,3.07836e-09,-2.21227e-13,12235,-14.7555], Tmin=(904.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC(F)=CC1(F)[CH][C]=C1(846)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,D} {10,S}
6  C u1 p0 c0 {4,S} {9,S} {12,S}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (0.97231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2750,3150,900,1100,182,240,577,636,1210,1413,180,680.118,680.131,680.138,680.143,680.144,680.161,680.194],'cm^-1')),
        HinderedRotor(inertia=(0.00851656,'amu*angstrom^2'), symmetry=1, barrier=(2.79562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46613,0.0584645,-5.06418e-05,2.29006e-08,-4.3419e-12,205.749,21.9351], Tmin=(100,'K'), Tmax=(1216.38,'K')), NASAPolynomial(coeffs=[10.5476,0.0286004,-1.38142e-05,2.71627e-09,-1.93429e-13,-2003.55,-23.6568], Tmin=(1216.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.97231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
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
    E0 = (40.8278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (211.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (74.3139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (214.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-84.8763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (121.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-13.1626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-91.1523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (137.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (352.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (295.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (401.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (101.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (127.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (353.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (320.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (483.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (180.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (77.3394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (136.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (117.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (157.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (85.8242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (66.6597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (162.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (111.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['C2HF(58)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(133.235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 133.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]C=C(F)C1C=C1F(828)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.6691e+10,'s^-1'), n=0.663, Ea=(162.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-OneDe;C]
Euclidian distance = 0
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
    reactants = ['F[C]=CC1=CC(F)(F)[CH]1(829)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(310.543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['FC1=CC2(F)C=C(F)C12(603)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;Cdsinglepri_rad_out]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['F[C]=CC1(F)C2[C](F)C21(830)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['FC1=CC2(F)[CH]C1(F)[CH]2(831)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.18, Ea=(79.245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_cdsingle]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['F[C]1[CH]C2(F)C=C(F)C12(832)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.2e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra_cdsingle] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_cdsingle]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(F)C=C(F)C=[C]F(833)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'F[C]=CC1=CC(F)=C1(834)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(13.7527,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]F(252)', 'FC1=CC(F)=C1(307)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[C]=CC1(F)C=C=C1(835)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(8.44811,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C2HF(58)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0816,'m^3/(mol*s)'), n=2.24, Ea=(3.06412,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'FC#CC1(F)[CH]C(F)=C1(836)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.29124e+14,'m^3/(mol*s)'), n=-1.78703, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]F(252)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', 'F[C]=CC1=C=C(F)[CH]1(837)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(153.349,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'F[C]=CC1(F)[C]=C=C1(838)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(206.5,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'F[C]=C=C1[CH]C(F)=C1(839)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(257.346,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['FC=[C]C1(F)[CH]C(F)=C1(840)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C]=CC1(F)[C]=C(F)C1(841)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]=[C]C1(F)C=C(F)C1(842)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3367.88,'s^-1'), n=2.74476, Ea=(136.665,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_23cy4;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC=CC1(F)[C]=C(F)[CH]1(843)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.07654e+10,'s^-1'), n=1.20849, Ea=(247.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_Cd;XH_out] + [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_single]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['F[C]C=C1C=C(F)C1F(844)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(178.232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    products = ['F[C](F)C=C1[CH]C(F)=C1(639)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(159.067,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C]=CC1(F)C=[C]C1F(845)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(132.695,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC(F)=CC1(F)[CH][C]=C1(846)'],
    products = ['F[C]=CC1(F)[CH]C(F)=C1(600)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(193.596,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #138',
    isomers = [
        'F[C]=CC1(F)[CH]C(F)=C1(600)',
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
    label = 'PDepNetwork #138',
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

