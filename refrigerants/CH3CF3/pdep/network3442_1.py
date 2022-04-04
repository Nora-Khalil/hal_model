species(
    label = 'F[C]=CC1[C](F)C1(F)F(11616)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-203.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,3010,987.5,1337.5,450,1655,167,640,1190,180.121,562.599,572.438,1486.03,2171.41],'cm^-1')),
        HinderedRotor(inertia=(0.00918512,'amu*angstrom^2'), symmetry=1, barrier=(2.16195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23831,0.0629503,-6.73969e-05,3.73946e-08,-8.41597e-12,-24377.3,27.559], Tmin=(100,'K'), Tmax=(1065.17,'K')), NASAPolynomial(coeffs=[11.8687,0.023031,-1.11823e-05,2.21152e-09,-1.58462e-13,-26641.9,-24.3976], Tmin=(1065.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC1=CC1(F)F(974)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,D} {7,S}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-329.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,323,467,575,827,1418,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2932.45,'J/mol'), sigma=(5.0155,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.04 K, Pc=52.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90533,0.00560985,9.08669e-05,-1.89954e-07,1.15706e-10,-39564.3,10.133], Tmin=(10,'K'), Tmax=(565.045,'K')), NASAPolynomial(coeffs=[4.14369,0.0253568,-1.84551e-05,6.16339e-09,-7.67829e-13,-39933.5,6.09131], Tmin=(565.045,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-329.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]=CC1(F)[CH]C1(F)F(11615)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {10,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-186.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,1230.04],'cm^-1')),
        HinderedRotor(inertia=(0.163159,'amu*angstrom^2'), symmetry=1, barrier=(3.75134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.991783,0.0723079,-0.000104572,8.56658e-08,-2.88866e-11,-22309.5,27.1478], Tmin=(100,'K'), Tmax=(720.402,'K')), NASAPolynomial(coeffs=[8.71306,0.0294352,-1.53023e-05,3.0539e-09,-2.17562e-13,-23422,-7.57062], Tmin=(720.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1=CC2C(F)(F)C12F(11654)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u0 p0 c0 {4,S} {6,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-512.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32079,0.0614859,-6.0196e-05,2.94883e-08,-5.88391e-12,-61602,18.8556], Tmin=(100,'K'), Tmax=(1183.26,'K')), NASAPolynomial(coeffs=[12.4994,0.0236967,-1.22908e-05,2.49767e-09,-1.8127e-13,-64247.4,-36.9555], Tmin=(1183.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCCCF) + group(CsCsCsFF) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + polycyclic(s2_3_4_ene_1)"""),
)

species(
    label = 'FC=CC1=C(F)C1(F)F(11655)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  C u0 p0 c0 {6,S} {9,D} {10,S}
9  C u0 p0 c0 {4,S} {8,D} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-481.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.4396,0.0703486,-7.50362e-05,3.83033e-08,-7.52981e-12,-57745.8,24.2878], Tmin=(100,'K'), Tmax=(1252.9,'K')), NASAPolynomial(coeffs=[18.8422,0.0115966,-4.69693e-06,8.75899e-10,-6.16416e-14,-62357.1,-68.643], Tmin=(1252.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-481.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F)) + ring(Cd-Cd-Cs(F)(F))"""),
)

species(
    label = 'F[C]=C[CH]C(F)=C(F)F(3442)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u1 p0 c0 {6,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-313.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,271,519,563,612,1379,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,167,640,1190,571.676,2933.91],'cm^-1')),
        HinderedRotor(inertia=(0.38152,'amu*angstrom^2'), symmetry=1, barrier=(88.4521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.381503,'amu*angstrom^2'), symmetry=1, barrier=(88.4435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.878125,0.0673906,-7.45023e-05,4.1691e-08,-9.22089e-12,-37549,26.6251], Tmin=(100,'K'), Tmax=(1101.37,'K')), NASAPolynomial(coeffs=[14.2863,0.0186942,-8.18055e-06,1.54597e-09,-1.08338e-13,-40502.4,-39.356], Tmin=(1101.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + group(CdCFH) + radical(C=CCJC=C) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'F[C]=CC1=C(F)C1(F)F(11656)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  C u0 p0 c0 {6,S} {9,D} {10,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-225.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5237,'amu*angstrom^2'), symmetry=1, barrier=(35.0329,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.682437,0.0739946,-9.56017e-05,6.0178e-08,-1.47487e-11,-27007.1,23.3489], Tmin=(100,'K'), Tmax=(1002.86,'K')), NASAPolynomial(coeffs=[15.4975,0.0149035,-7.21784e-06,1.42353e-09,-1.01947e-13,-29978.6,-48.1674], Tmin=(1002.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cd-Cs(F)(F)) + radical(Cdj(Cd-CdH)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F))"""),
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
    label = 'F[C]=CC1C(F)=C1F(11657)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {4,S} {5,D}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (75.5671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,260,386,409,525,515,635,761,893,1354,1482,3010,987.5,1337.5,450,1655,167,640,1190,180,180,2101.12],'cm^-1')),
        HinderedRotor(inertia=(0.00383931,'amu*angstrom^2'), symmetry=1, barrier=(11.8423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27225,0.0647913,-9.65056e-05,7.99318e-08,-2.66581e-11,9182.27,21.6656], Tmin=(100,'K'), Tmax=(789.226,'K')), NASAPolynomial(coeffs=[8.62345,0.0238191,-1.15744e-05,2.22621e-09,-1.54593e-13,8137.6,-11.3268], Tmin=(789.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.5671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cd(F)-Cd) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'F[C]1[CH]C1(F)F(11362)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-103.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2950,1000,212,367,445,1450,2726.48,2726.56],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00781,0.0264579,-1.29283e-05,-2.29098e-08,2.96552e-11,-12413.9,18.9514], Tmin=(100,'K'), Tmax=(451.813,'K')), NASAPolynomial(coeffs=[4.42599,0.0197907,-1.03421e-05,2.11919e-09,-1.54582e-13,-12602.1,12.5711], Tmin=(451.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FC#CC1[C](F)C1(F)F(11658)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,T}
9  C u0 p0 c0 {4,S} {8,T}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-196.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,2175,525,239,401,1367,180,1311.54,1311.55,1311.56,1311.6],'cm^-1')),
        HinderedRotor(inertia=(0.0202603,'amu*angstrom^2'), symmetry=1, barrier=(5.02936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35278,0.0623302,-8.02661e-05,5.64609e-08,-1.61732e-11,-23558.5,25.9374], Tmin=(100,'K'), Tmax=(846.601,'K')), NASAPolynomial(coeffs=[9.58237,0.0234495,-1.13815e-05,2.22007e-09,-1.56879e-13,-24952,-12.3957], Tmin=(846.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Ct-CtCs) + group(CtCF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'F[C]=CC1=C(F)[C]1F(11659)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {8,S}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u1 p0 c0 {2,S} {4,S} {6,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 C u0 p0 c0 {4,S} {8,D} {9,S}
8 C u1 p0 c0 {3,S} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (201.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,271,519,563,612,1379,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24577,'amu*angstrom^2'), symmetry=1, barrier=(28.6426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04176,0.0653607,-8.54948e-05,5.41838e-08,-1.3301e-11,24375.8,22.1276], Tmin=(100,'K'), Tmax=(1004.16,'K')), NASAPolynomial(coeffs=[14.4881,0.0117976,-5.48175e-06,1.06198e-09,-7.53758e-14,21675.4,-42.7987], Tmin=(1004.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + ring(Cd-Cs(F)-Cd) + radical(CsCdCdF1s_ring) + radical(Cdj(Cd-CdH)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C]=C[C]1C(F)C1(F)F(3619)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {11,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-271.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40137,0.05856,-5.56587e-05,2.75591e-08,-5.60886e-12,-32505.6,23.172], Tmin=(100,'K'), Tmax=(1160.22,'K')), NASAPolynomial(coeffs=[11.2707,0.0245345,-1.16687e-05,2.28241e-09,-1.62365e-13,-34795.7,-25.9086], Tmin=(1160.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(Allyl_T) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FC=[C]C1[C](F)C1(F)F(11660)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-215.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25685,0.0642257,-7.31779e-05,4.4422e-08,-1.10395e-11,-25855.9,27.481], Tmin=(100,'K'), Tmax=(965.706,'K')), NASAPolynomial(coeffs=[10.7034,0.0250963,-1.23974e-05,2.46136e-09,-1.7647e-13,-27680.3,-17.7633], Tmin=(965.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Cds_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=[C]C1C(F)C1(F)F(11661)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {5,S} {6,S} {11,S}
8  C u1 p0 c0 {5,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-156.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,1685,370,167,640,1190,180.069,938.013,938.786,938.87,2566.84],'cm^-1')),
        HinderedRotor(inertia=(0.232671,'amu*angstrom^2'), symmetry=1, barrier=(20.3789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0421,0.0649569,-6.91858e-05,3.71278e-08,-7.94335e-12,-18671.3,26.0294], Tmin=(100,'K'), Tmax=(1129.03,'K')), NASAPolynomial(coeffs=[13.7339,0.0199915,-9.44576e-06,1.85249e-09,-1.3236e-13,-21537.1,-36.7413], Tmin=(1129.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C[C]1[C](F)C1(F)F(11649)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {7,S} {8,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {6,S} {9,D} {10,S}
9  C u0 p0 c0 {4,S} {8,D} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-321.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,212,367,445,1450,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,1022.67,1022.82,1026.63],'cm^-1')),
        HinderedRotor(inertia=(0.00496134,'amu*angstrom^2'), symmetry=1, barrier=(3.69165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38818,0.0565035,-5.14908e-05,2.41096e-08,-4.5894e-12,-38588.8,25.4698], Tmin=(100,'K'), Tmax=(1244.77,'K')), NASAPolynomial(coeffs=[12.0786,0.0221503,-1.00936e-05,1.93821e-09,-1.3647e-13,-41250.2,-28.4459], Tmin=(1244.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(Allyl_T) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]1[C](F)C1C=C(F)F(11662)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {7,S}
7  C u1 p0 c0 {2,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-253.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60945,0.0597188,-6.59189e-05,3.2652e-08,-7.00233e-14,-30402.8,29.3238], Tmin=(100,'K'), Tmax=(560.093,'K')), NASAPolynomial(coeffs=[6.82875,0.0317909,-1.6156e-05,3.21487e-09,-2.29555e-13,-31134.1,5.8603], Tmin=(560.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-112.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (62.3198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-103.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-48.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (9.41822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (114.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (79.7472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (264.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (85.9944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (106.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (338.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (180.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (45.6162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (57.5968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (67.8505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (73.1396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (105.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['C2HF(58)', 'FC1=CC1(F)F(974)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]=CC1(F)[CH]C1(F)F(11615)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['FC1=CC2C(F)(F)C12F(11654)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_noH;Cdsinglepri_rad_out]
Euclidian distance = 1.7320508075688772
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['FC=CC1=C(F)C1(F)F(11655)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]=C[CH]C(F)=C(F)F(3442)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]F(489)', 'FC1=CC1(F)F(974)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(1.54577,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'F[C]=CC1=C(F)C1(F)F(11656)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(2.1335,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'F[C]=CC1C(F)=C1F(11657)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(24.3867,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2HF(58)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0408,'m^3/(mol*s)'), n=2.24, Ea=(2.78971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_N-Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'FC#CC1[C](F)C1(F)F(11658)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.29124e+14,'m^3/(mol*s)'), n=-1.78703, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]F(489)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['HF(38)', 'F[C]=CC1=C(F)[C]1F(11659)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.28222,'m^3/(mol*s)'), n=1.29695, Ea=(168.681,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['F[C]=C[C]1C(F)C1(F)F(3619)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.71585e+08,'s^-1'), n=1.23602, Ea=(157.766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_noH;Cs_H_out_noH] + [R2H_S_cy3;C_rad_out_single;Cs_H_out_noH] + [R2H_S_cy3;C_rad_out_noH;Cs_H_out] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['FC=[C]C1[C](F)C1(F)F(11660)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[C]=[C]C1C(F)C1(F)F(11661)'],
    products = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.1e+10,'s^-1'), n=0.78, Ea=(132.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_23cy3;Cd_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['FC=C[C]1[C](F)C1(F)F(11649)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(323312,'s^-1'), n=2.11016, Ea=(185.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C]=CC1[C](F)C1(F)F(11616)'],
    products = ['F[C]1[C](F)C1C=C(F)F(11662)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(217.623,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #3442',
    isomers = [
        'F[C]=CC1[C](F)C1(F)F(11616)',
    ],
    reactants = [
        ('C2HF(58)', 'FC1=CC1(F)F(974)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3442',
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

