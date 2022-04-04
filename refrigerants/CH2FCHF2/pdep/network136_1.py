species(
    label = '[CH]=C(F)C1(F)[CH]C(F)=C1(598)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {6,D}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-17.4277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,246,474,533,1155,3120,650,792.5,1650,623.479,624.331,625.298,625.353,625.448,625.856],'cm^-1')),
        HinderedRotor(inertia=(0.00918332,'amu*angstrom^2'), symmetry=1, barrier=(2.55005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15382,0.0649699,-6.35072e-05,3.23164e-08,-6.75995e-12,-1995.56,22.2378], Tmin=(100,'K'), Tmax=(1130.33,'K')), NASAPolynomial(coeffs=[11.8864,0.0269898,-1.31064e-05,2.59033e-09,-1.85373e-13,-4421.86,-30.8557], Tmin=(1130.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.4277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    label = '[CH]C(F)=C(F)C1C=C1F(2481)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {9,S}
9  C u2 p0 c0 {8,S} {12,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (124.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,260,323,386,409,467,525,515,575,635,761,827,893,1354,1418,1482,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839809,0.0772716,-8.9575e-05,5.11378e-08,-6.87894e-12,15110.8,24.8537], Tmin=(100,'K'), Tmax=(601.68,'K')), NASAPolynomial(coeffs=[8.61314,0.0370683,-1.79523e-05,3.474e-09,-2.4372e-13,13967.7,-10.425], Tmin=(601.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs-Cd(F)-Cd) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(F)C1[C](F)C=C1F(599)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-12.9964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,246,474,533,1155,3120,650,792.5,1650,180,180,180,180,1054.92,1054.92,1054.92,1055.02],'cm^-1')),
        HinderedRotor(inertia=(0.00540575,'amu*angstrom^2'), symmetry=1, barrier=(4.26925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3431.89,'J/mol'), sigma=(5.38731,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.05 K, Pc=49.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936382,0.0678853,-6.81855e-05,3.50795e-08,-7.29323e-12,-1453.19,24.0022], Tmin=(100,'K'), Tmax=(1150.4,'K')), NASAPolynomial(coeffs=[13.4279,0.0244515,-1.15523e-05,2.25995e-09,-1.61016e-13,-4327.23,-38.0122], Tmin=(1150.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.9964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=C(F)C1=CC(F)(F)[CH]1(2482)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
5  C u0 p0 c0 {6,S} {7,D} {8,S}
6  C u1 p0 c0 {4,S} {5,S} {11,S}
7  C u0 p0 c0 {4,S} {5,D} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-8.74064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2750,3150,900,1100,250,446,589,854,899,3120,650,792.5,1650,180,815.518,816.028,816.757,816.869,816.951,816.975,817.393],'cm^-1')),
        HinderedRotor(inertia=(0.109383,'amu*angstrom^2'), symmetry=1, barrier=(2.51494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615851,0.0639532,-4.87952e-05,1.24016e-08,6.99486e-13,-920.466,22.5611], Tmin=(100,'K'), Tmax=(1085.71,'K')), NASAPolynomial(coeffs=[17.9651,0.017047,-7.49398e-06,1.4732e-09,-1.07253e-13,-5690.4,-67.1833], Tmin=(1085.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.74064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFF) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-CdsHH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(cyclobutene-allyl) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[CH]=C(F)C1(F)C2[C](F)C21(2483)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u1 p0 c0 {2,S} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (65.2031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969924,0.0736076,-9.52958e-05,6.98066e-08,-2.13963e-11,7945,23.282], Tmin=(100,'K'), Tmax=(783.866,'K')), NASAPolynomial(coeffs=[8.9507,0.0328808,-1.73581e-05,3.51894e-09,-2.54172e-13,6693.88,-13.2769], Tmin=(783.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.2031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFH) + group(CdCsCdF) + group(Cds-CdsHH) + polycyclic(s2_3_3_ane) + radical(CsCsCsF1s) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    label = '[CH]=C(F)C=C(F)C(=[CH])F(2484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,D} {7,S}
5  C u0 p0 c0 {4,D} {6,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {6,D} {11,S}
9  C u1 p0 c0 {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (92.6229,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,3010,987.5,1337.5,450,1655,174,326,386,506,510,668,789,919,801,997,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(0.857396,'amu*angstrom^2'), symmetry=1, barrier=(19.7132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857091,'amu*angstrom^2'), symmetry=1, barrier=(19.7062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271093,0.0831151,-0.000101067,6.06432e-08,-1.42899e-11,11273.4,26.725], Tmin=(100,'K'), Tmax=(1037.32,'K')), NASAPolynomial(coeffs=[16.6564,0.0199319,-9.70252e-06,1.9246e-09,-1.3838e-13,7874.09,-52.925], Tmin=(1037.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.6229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)-Cds(F)) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)-Cds(F)) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cdj(Cd-CdF1s)(H)) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[CH]=C(F)C1=CC(F)=C1(2485)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,D} {5,S} {7,S}
4  C u0 p0 c0 {3,D} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u0 p0 c0 {2,S} {3,S} {8,D}
8  C u1 p0 c0 {7,D} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (353.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,280,518,736,852,873,250,446,589,854,899,3120,650,792.5,1650,180,180,1274.65,1274.72,1274.72,1274.73,1274.75,1274.9],'cm^-1')),
        HinderedRotor(inertia=(0.285193,'amu*angstrom^2'), symmetry=1, barrier=(6.55715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.547,0.0526692,-4.4655e-05,1.92158e-08,-3.37135e-12,42643.1,20.8289], Tmin=(100,'K'), Tmax=(1335.95,'K')), NASAPolynomial(coeffs=[11.8567,0.0218012,-9.99699e-06,1.92096e-09,-1.34972e-13,39888.4,-31.8957], Tmin=(1335.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(CdCCF) + group(Cds-CdsHH) + ring(Cd-Cd-Cd-Cd(F)) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[CH]=C(F)C1(F)C=C=C1(2486)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u0 p0 c0 {4,D} {5,D}
8  C u1 p0 c0 {6,D} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (395.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,246,474,533,1155,3120,650,792.5,1650,180,180,180,180,1386.49,1388.35,1391.92,1393.4],'cm^-1')),
        HinderedRotor(inertia=(0.353259,'amu*angstrom^2'), symmetry=1, barrier=(8.12211,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00928,0.0765562,-0.000134529,1.30787e-07,-4.9079e-11,47612.1,21.5026], Tmin=(100,'K'), Tmax=(820.653,'K')), NASAPolynomial(coeffs=[4.78769,0.0357644,-1.90709e-05,3.7703e-09,-2.63963e-13,47745.4,8.61132], Tmin=(820.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'C#CC1(F)[CH]C(F)=C1(2487)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {9,S}
6  C u0 p0 c0 {2,S} {4,S} {5,D}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (98.7602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,2175,525,750,770,3400,2100,386.328,386.569,386.744,386.794,386.813],'cm^-1')),
        HinderedRotor(inertia=(0.00112954,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22351,0.0563997,-5.19675e-05,2.4706e-08,-4.69039e-12,11982,18.8778], Tmin=(100,'K'), Tmax=(1269.87,'K')), NASAPolynomial(coeffs=[13.1821,0.0187313,-7.47295e-06,1.34703e-09,-9.1732e-14,8944.77,-41.6726], Tmin=(1269.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.7602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl)"""),
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
    label = 'F2(78)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C=C1[CH]C(F)=C1(2488)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {6,D}
3  C u1 p0 c0 {2,S} {5,S} {9,S}
4  C u0 p0 c0 {2,S} {5,D} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {4,D}
6  C u0 p0 c0 {2,D} {7,D}
7  C u1 p0 c0 {6,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (408.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,271,519,563,612,1379,540,610,2055,3120,650,792.5,1650,652.491,654.565,654.791,655.129,655.813,656.111,657.349,659.999],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0862,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8261,0.0279896,4.52102e-05,-9.5018e-08,4.35924e-11,49168.6,16.8782], Tmin=(100,'K'), Tmax=(906.569,'K')), NASAPolynomial(coeffs=[19.9043,0.000602571,3.85986e-06,-8.79301e-10,5.75641e-14,43738.4,-80.4364], Tmin=(906.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(C=C=CJ)"""),
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
    label = '[CH]=C(F)C1=C=C(F)[CH]1(2489)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {6,S} {7,D}
4  C u1 p0 c0 {3,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u0 p0 c0 {3,D} {5,D}
8  C u1 p0 c0 {6,D} {10,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (536.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,132,421,945,1169,1268,250,446,589,854,899,3120,650,792.5,1650,614.899,614.983,615.154,615.336,615.379,615.579,615.621],'cm^-1')),
        HinderedRotor(inertia=(0.0122879,'amu*angstrom^2'), symmetry=1, barrier=(3.30633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24515,0.0552556,-5.20697e-05,2.41442e-08,-4.40966e-12,64591.6,22.0362], Tmin=(100,'K'), Tmax=(1324.59,'K')), NASAPolynomial(coeffs=[14.5525,0.0150699,-6.56209e-06,1.24014e-09,-8.67691e-14,61066.2,-45.9046], Tmin=(1324.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCddCF) + group(CdCCF) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(C=CCJC=C) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = '[CH]=C(F)C1(F)[C]=C=C1(2490)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {8,D}
6  C u1 p0 c0 {3,S} {7,D}
7  C u0 p0 c0 {4,D} {6,D}
8  C u1 p0 c0 {5,D} {10,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (632.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,246,474,533,1155,3120,650,792.5,1650,180,180,180,180,1526.62,1527.15,1528.49],'cm^-1')),
        HinderedRotor(inertia=(0.400595,'amu*angstrom^2'), symmetry=1, barrier=(9.21047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.956795,0.0813789,-0.000160028,1.60786e-07,-6.02145e-11,76216.3,22.8579], Tmin=(100,'K'), Tmax=(847.339,'K')), NASAPolynomial(coeffs=[3.75929,0.034983,-1.91831e-05,3.77966e-09,-2.61966e-13,76932,16.8276], Tmin=(847.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cds_S) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=C(F)C1(F)[C]=C(F)C1(2491)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {4,S} {6,D}
9  C u1 p0 c0 {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (70.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,161,331,427,521,441,625,1076,1234,3120,650,792.5,1650,320.632,320.673,320.768,320.776,1427.61,1990.98,1991.14],'cm^-1')),
        HinderedRotor(inertia=(0.418631,'amu*angstrom^2'), symmetry=1, barrier=(30.5438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816056,0.0759201,-0.00010069,7.41311e-08,-2.23932e-11,8590.34,25.0399], Tmin=(100,'K'), Tmax=(801.887,'K')), NASAPolynomial(coeffs=[9.98475,0.0301815,-1.51263e-05,2.99123e-09,-2.12803e-13,7119.99,-17.1688], Tmin=(801.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-vinyl) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'C=C(F)C1(F)[C]=C(F)[CH]1(2492)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u1 p0 c0 {4,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {3,S} {5,S} {9,D}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-23.2026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19689,0.0664845,-7.1304e-05,4.17903e-08,-1.02379e-11,-2693.71,21.943], Tmin=(100,'K'), Tmax=(968.089,'K')), NASAPolynomial(coeffs=[10.002,0.0301034,-1.49336e-05,2.97128e-09,-2.13281e-13,-4398.53,-20.2509], Tmin=(968.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]C(F)=C1C=C(F)C1F(2493)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {8,D}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {3,S} {5,D} {9,S}
9  C u2 p0 c0 {8,S} {12,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-50.0175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.319223,0.0718057,-6.32713e-05,2.76419e-08,-4.77795e-12,-5875.6,24.5158], Tmin=(100,'K'), Tmax=(1393.64,'K')), NASAPolynomial(coeffs=[17.9229,0.0212802,-8.89004e-06,1.628e-09,-1.11437e-13,-10782.3,-66.2549], Tmin=(1393.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.0175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(AllylJ2_triplet) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[CH]C(F)=C1[CH]C(F)=C1(682)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {8,D}
5  C u1 p0 c0 {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,D}
8  C u0 p0 c0 {2,S} {4,D} {9,S}
9  C u1 p0 c0 {3,S} {8,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-152.453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,206,336,431,607,515,611,528,696,1312,1446,234,589,736,816,1240,3237,693.968,694.626,694.878,695.452,695.509,695.511,695.614,695.651,696.575],'cm^-1')),
        HinderedRotor(inertia=(0.000348127,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811884,0.0521149,1.33193e-06,-5.29468e-08,2.77661e-11,-18204.6,22.3293], Tmin=(100,'K'), Tmax=(939.739,'K')), NASAPolynomial(coeffs=[21.7109,0.00845796,-1.28993e-06,2.08642e-10,-2.11149e-14,-24132.7,-87.8397], Tmin=(939.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = '[CH]=C(F)C1(F)C=[C]C1F(2494)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {5,S} {6,D}
9  C u1 p0 c0 {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (104.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,2950,1000,246,474,533,1155,3120,650,792.5,1650,180,180,180,959.779,963.578],'cm^-1')),
        HinderedRotor(inertia=(0.00390132,'amu*angstrom^2'), symmetry=1, barrier=(2.57246,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05887,0.0702823,-8.2098e-05,5.25856e-08,-1.39879e-11,12716.4,24.4229], Tmin=(100,'K'), Tmax=(899.279,'K')), NASAPolynomial(coeffs=[10.0634,0.0302307,-1.52926e-05,3.06105e-09,-2.20192e-13,11096.8,-18.0629], Tmin=(899.279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + radical(Cdj(Cd-CsF1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C(F)C1(F)[CH][C]=C1(2495)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {8,D}
6  C u1 p0 c0 {4,S} {9,S} {11,S}
7  C u0 p0 c0 {4,S} {9,D} {10,S}
8  C u0 p0 c0 {3,S} {5,D} {12,S}
9  C u1 p0 c0 {6,S} {7,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (18.4252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,323,467,575,827,1418,2750,3150,900,1100,194,682,905,1196,1383,3221,180,180,180,917.64,917.758,917.764,917.852,917.907],'cm^-1')),
        HinderedRotor(inertia=(0.091944,'amu*angstrom^2'), symmetry=1, barrier=(2.11397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.272,0.0652246,-6.92016e-05,4.07015e-08,-1.01073e-11,2309.91,21.6551], Tmin=(100,'K'), Tmax=(950.912,'K')), NASAPolynomial(coeffs=[9.35434,0.031227,-1.55737e-05,3.10466e-09,-2.23052e-13,772.763,-16.9309], Tmin=(950.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.4252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]=[C]C1(F)C=C(F)C1F(2496)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  C u1 p0 c0 {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (93.2533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,2950,1000,323,467,575,827,1418,1685,370,3120,650,792.5,1650,180,180,1646.91],'cm^-1')),
        HinderedRotor(inertia=(0.860398,'amu*angstrom^2'), symmetry=1, barrier=(19.7822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987166,0.072235,-9.20734e-05,6.61685e-08,-1.97423e-11,11318.9,25.7841], Tmin=(100,'K'), Tmax=(807.469,'K')), NASAPolynomial(coeffs=[9.28355,0.0311353,-1.57212e-05,3.12782e-09,-2.23598e-13,9979.18,-12.4667], Tmin=(807.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.2533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Cds_S) + radical(Cds_P) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
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
    E0 = (32.0351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (194.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (56.9562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (203.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-102.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (103.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-30.5202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-93.6254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (138.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (346.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (287.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (384.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (141.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (92.2276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (344.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (324.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (311.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (465.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (119.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-65.4566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (56.9568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (55.3122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (145.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (106.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (176.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (76.5653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['C2HF(58)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(141.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 141.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(F)=C(F)C1C=C1F(2481)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.6691e+10,'s^-1'), n=0.663, Ea=(162.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-OneDe;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(F)C1[C](F)C=C1F(599)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.13382e+11,'s^-1'), n=0.663, Ea=(162.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ-OneDe;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(F)C1=CC(F)(F)[CH]1(2482)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(304.474,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['FC1=CC2(F)C(F)=CC12(602)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['[CH]=C(F)C1(F)C2[C](F)C21(2483)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['FC1=CC2(F)[CH]C1(F)[CH]2(831)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.18, Ea=(79.245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['F[C]1[CH]C2(F)C(F)=CC12(814)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.12675e+11,'s^-1'), n=0.175, Ea=(16.1398,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_cyclic;multiplebond_intra;radadd_intra_cdsingleH] + [R6_cyclic;doublebond_intra;radadd_intra_cdsingle] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(F)C=C(F)C(=[CH])F(2484)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH]=C(F)C1=CC(F)=C1(2485)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(12.4517,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]F(252)', 'FC1=CC(F)=C1(307)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.0102528,'m^3/(mol*s)'), n=2.42326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7837337452902932, var=8.455757054931007, Tref=1000.0, N=46, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[CH]=C(F)C1(F)C=C=C1(2486)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(8.44811,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C#CC1(F)[CH]C(F)=C1(2487)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.0689,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C2HF(58)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0816,'m^3/(mol*s)'), n=2.24, Ea=(2.42856,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]F(252)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(78)', '[CH]=C=C1[CH]C(F)=C1(2488)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(17.4436,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', '[CH]=C(F)C1=C=C(F)[CH]1(2489)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(149.071,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', '[CH]=C(F)C1(F)[C]=C=C1(2490)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(206.5,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(F)C1(F)[C]=C(F)C1(2491)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['C=C(F)C1(F)[C]=C(F)[CH]1(2492)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['[CH]C(F)=C1C=C(F)C1F(2493)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(166.722,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['F[CH]C(F)=C1[CH]C(F)=C1(682)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(165.077,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(F)C1(F)C=[C]C1F(2494)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(132.695,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FC=C(F)C1(F)[CH][C]=C1(2495)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(180.752,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]C1(F)C=C(F)C1F(2496)'],
    products = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(175.841,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(F)C1(F)[CH]C(F)=C1(598)'],
    products = ['FC=[C]C1(F)[CH]C(F)=C1(840)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(186.33,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #136',
    isomers = [
        '[CH]=C(F)C1(F)[CH]C(F)=C1(598)',
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
    label = 'PDepNetwork #136',
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

