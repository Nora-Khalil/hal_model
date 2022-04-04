species(
    label = '[CH]=C(F)C(F)[CH]C(F)(F)F(12108)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u1 p0 c0 {9,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-643.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,246,474,533,1155,3120,650,792.5,1650,180,1593.9],'cm^-1')),
        HinderedRotor(inertia=(0.202529,'amu*angstrom^2'), symmetry=1, barrier=(4.65655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203387,'amu*angstrom^2'), symmetry=1, barrier=(4.67627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20128,'amu*angstrom^2'), symmetry=1, barrier=(4.62781,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0690735,0.0949013,-0.000150342,1.31195e-07,-4.59972e-11,-77213.8,32.0729], Tmin=(100,'K'), Tmax=(771.117,'K')), NASAPolynomial(coeffs=[10.0222,0.0352077,-1.85377e-05,3.68176e-09,-2.60065e-13,-78509.1,-11.8035], Tmin=(771.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFF) + group(CdCsCdF) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    label = 'FC=CC(F)(F)F(1027)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {8,S}
7 C u0 p0 c0 {4,S} {6,D} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-831.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,224.569],'cm^-1')),
        HinderedRotor(inertia=(0.0840844,'amu*angstrom^2'), symmetry=1, barrier=(1.93327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3084.67,'J/mol'), sigma=(4.9,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77934,0.0173102,9.62417e-05,-2.61573e-07,1.94701e-10,-100018,11.5092], Tmin=(10,'K'), Tmax=(474.505,'K')), NASAPolynomial(coeffs=[4.81364,0.0316651,-2.20773e-05,7.14072e-09,-8.67875e-13,-100376,4.55317], Tmin=(474.505,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-831.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FCDCC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=C(F)C([CH]F)C(F)(F)F(6905)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u1 p0 c0 {5,S} {6,S} {12,S}
10 C u1 p0 c0 {8,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-643.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,246,474,533,1155,334,575,1197,1424,3202,3120,650,792.5,1650,279.03,283.678,2833.83],'cm^-1')),
        HinderedRotor(inertia=(0.00208074,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247797,'amu*angstrom^2'), symmetry=1, barrier=(14.3859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638935,'amu*angstrom^2'), symmetry=1, barrier=(36.465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3176.14,'J/mol'), sigma=(5.55441,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=496.11 K, Pc=42.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0761416,0.0916738,-0.000122794,8.50211e-08,-2.34755e-11,-77226.5,31.7243], Tmin=(100,'K'), Tmax=(884.01,'K')), NASAPolynomial(coeffs=[14.2932,0.0273414,-1.36288e-05,2.69215e-09,-1.91683e-13,-79740,-35.1113], Tmin=(884.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCsFHH) + group(CdCsCdF) + group(Cds-CdsHH) + radical(CsCsF1sH) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]C(F)(F)F(462)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 C u2 p0 c0 {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-319.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([193,295,551,588,656,1146,1192,1350,180,747.919,1907.62],'cm^-1')),
        HinderedRotor(inertia=(0.00850025,'amu*angstrom^2'), symmetry=1, barrier=(3.32661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44977,0.0297269,-2.96856e-05,1.40699e-08,-2.55624e-12,-38337.5,14.102], Tmin=(100,'K'), Tmax=(1355.02,'K')), NASAPolynomial(coeffs=[10.9803,0.0045449,-1.80947e-06,3.54966e-10,-2.58626e-14,-40649.3,-29.6449], Tmin=(1355.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(F)=CF(878)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 C u2 p0 c0 {3,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-0.651196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([323,467,575,827,1418,194,682,905,1196,1383,3221,439.435,439.564,439.597],'cm^-1')),
        HinderedRotor(inertia=(0.377592,'amu*angstrom^2'), symmetry=1, barrier=(51.7341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59677,0.0297955,-1.83591e-05,5.50976e-09,-6.73981e-13,-27.3053,15.0855], Tmin=(100,'K'), Tmax=(1812.93,'K')), NASAPolynomial(coeffs=[9.00424,0.0156583,-6.66213e-06,1.20844e-09,-8.08368e-14,-2350.57,-19.6389], Tmin=(1812.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.651196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = 'FC1=CC(C(F)(F)F)C1F(12175)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
9  C u0 p0 c0 {6,S} {10,D} {13,S}
10 C u0 p0 c0 {5,S} {7,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-941.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679044,0.0655082,-5.31952e-05,2.06325e-08,-3.16424e-12,-113126,25.2337], Tmin=(100,'K'), Tmax=(1547,'K')), NASAPolynomial(coeffs=[18.0582,0.0205714,-9.62333e-06,1.85538e-09,-1.29774e-13,-118503,-66.1932], Tmin=(1547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-941.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCCFH) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(Cd(F)-Cd-Cs-Cs)"""),
)

species(
    label = 'C=C(F)C(F)=CC(F)(F)F(12209)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {4,S} {7,D} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,D}
10 C u0 p0 c0 {9,D} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-962.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446585,0.0757332,-7.49806e-05,3.67038e-08,-7.14203e-12,-115675,27.3229], Tmin=(100,'K'), Tmax=(1236.34,'K')), NASAPolynomial(coeffs=[16.7749,0.0229047,-1.08853e-05,2.14155e-09,-1.53139e-13,-119712,-54.9158], Tmin=(1236.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-962.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cds-CdsCsH) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)-Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)-Cds(F)) + group(Cds-CdsHH)"""),
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
    label = '[CH]=C(F)C=CC(F)(F)F(12210)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {6,D} {8,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-516.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,2995,3025,975,1000,1300,1375,400,500,1630,1680,250,446,589,854,899,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.586688,'amu*angstrom^2'), symmetry=1, barrier=(13.4891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.586677,'amu*angstrom^2'), symmetry=1, barrier=(13.4889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73775,0.0739117,-8.37371e-05,4.82168e-08,-1.11126e-11,-62029.6,24.7614], Tmin=(100,'K'), Tmax=(1048.95,'K')), NASAPolynomial(coeffs=[13.9949,0.0233594,-1.14495e-05,2.27545e-09,-1.63567e-13,-64810.9,-39.8306], Tmin=(1048.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-516.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-CdsHH) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[CH]=C(F)C(F)=CC(F)(F)F(12211)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {4,S} {7,D} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,D}
10 C u1 p0 c0 {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-699.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,280,518,736,852,873,250,446,589,854,899,3120,650,792.5,1650,302.074],'cm^-1')),
        HinderedRotor(inertia=(0.114106,'amu*angstrom^2'), symmetry=1, barrier=(2.62353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.557581,'amu*angstrom^2'), symmetry=1, barrier=(12.8199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442339,0.0794503,-9.12573e-05,5.20905e-08,-1.17791e-11,-84028.8,27.6466], Tmin=(100,'K'), Tmax=(1073.82,'K')), NASAPolynomial(coeffs=[15.896,0.0218852,-1.08455e-05,2.16783e-09,-1.5645e-13,-87347.6,-48.0088], Tmin=(1073.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cds-CdsCsH) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)-Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)-Cds(F)) + group(Cds-CdsHH) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = '[CH]=C(F)C(F)C=C(F)F(9342)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {4,S} {6,D}
9  C u1 p0 c0 {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-446.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,246,474,533,1155,182,240,577,636,1210,1413,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.469416,'amu*angstrom^2'), symmetry=1, barrier=(10.7928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39714,'amu*angstrom^2'), symmetry=1, barrier=(32.123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811748,0.0780303,-0.000104871,7.18106e-08,-1.67555e-11,-53596.8,26.2219], Tmin=(100,'K'), Tmax=(612.96,'K')), NASAPolynomial(coeffs=[9.62534,0.0310865,-1.58621e-05,3.13856e-09,-2.22505e-13,-54875.8,-13.6047], Tmin=(612.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-446.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCFF) + group(Cds-CdsHH) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'C#CC(F)[CH]C(F)(F)F(12212)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {6,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u0 p0 c0 {5,S} {9,T}
9  C u0 p0 c0 {8,T} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-537.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([323,546,651,945,1004,1329,3278,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.123371,'amu*angstrom^2'), symmetry=1, barrier=(2.83655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123249,'amu*angstrom^2'), symmetry=1, barrier=(2.83374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.0982,'amu*angstrom^2'), symmetry=1, barrier=(71.2337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504274,0.0834674,-0.00012838,1.08326e-07,-3.65591e-11,-64495.1,27.0634], Tmin=(100,'K'), Tmax=(795.713,'K')), NASAPolynomial(coeffs=[9.9142,0.029671,-1.47276e-05,2.84944e-09,-1.98186e-13,-65787.1,-14.892], Tmin=(795.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cs_S)"""),
)

species(
    label = 'F[CH][CH]C(F)(F)F(9334)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u1 p0 c0 {5,S} {7,S} {8,S}
7 C u1 p0 c0 {4,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-581.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,334,575,1197,1424,3202,208.405,1612.5],'cm^-1')),
        HinderedRotor(inertia=(0.236761,'amu*angstrom^2'), symmetry=1, barrier=(6.44327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00345971,'amu*angstrom^2'), symmetry=1, barrier=(6.42129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9911,0.0487013,-6.31668e-05,4.72978e-08,-1.48179e-11,-69825.2,23.0747], Tmin=(100,'K'), Tmax=(768.081,'K')), NASAPolynomial(coeffs=[7.09051,0.0221438,-1.13006e-05,2.27841e-09,-1.64263e-13,-70608.5,-0.181341], Tmin=(768.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-581.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    collisionModel = TransportData(shapeIndex=1, epsilon=(1045.13,'J/mol'), sigma=(3.301,'angstroms'), dipoleMoment=(0,'De'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C=C[CH]C(F)(F)F(12213)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5  C u1 p0 c0 {4,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {6,D} {8,D}
8  C u1 p0 c0 {7,D} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-237.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,311.538],'cm^-1')),
        HinderedRotor(inertia=(0.788612,'amu*angstrom^2'), symmetry=1, barrier=(18.1318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.8993,'amu*angstrom^2'), symmetry=1, barrier=(66.6607,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15338,0.0610427,-6.4068e-05,3.48142e-08,-7.53936e-12,-28478.9,23.0136], Tmin=(100,'K'), Tmax=(1120.14,'K')), NASAPolynomial(coeffs=[12.8595,0.0192406,-8.09027e-06,1.49845e-09,-1.03756e-13,-31101.4,-34.7897], Tmin=(1120.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(C=C=CJ)"""),
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
    label = '[CH]=C=C(F)[CH]C(F)(F)F(12214)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u1 p0 c0 {5,S} {7,S} {10,S}
7  C u0 p0 c0 {4,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {9,D}
9  C u1 p0 c0 {8,D} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-422.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,132,421,945,1169,1268,540,610,2055,3120,650,792.5,1650,417.546],'cm^-1')),
        HinderedRotor(inertia=(0.488749,'amu*angstrom^2'), symmetry=1, barrier=(11.2373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.68458,'amu*angstrom^2'), symmetry=1, barrier=(61.7238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.853095,0.0694745,-8.11947e-05,4.84904e-08,-1.1462e-11,-50695.1,25.3202], Tmin=(100,'K'), Tmax=(1032.32,'K')), NASAPolynomial(coeffs=[13.7395,0.0195417,-8.63886e-06,1.63319e-09,-1.14235e-13,-53355.7,-37.2587], Tmin=(1032.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-422.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CdCddCF) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(F)=C(F)CC(F)(F)F(12215)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {10,S}
10 C u2 p0 c0 {9,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-725.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266325,0.0875919,-0.000102288,6.54081e-08,-1.72439e-11,-87076.3,29.733], Tmin=(100,'K'), Tmax=(911.339,'K')), NASAPolynomial(coeffs=[11.9921,0.0361249,-1.75752e-05,3.43796e-09,-2.43916e-13,-89213.5,-25.7484], Tmin=(911.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-725.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(F)=C(F)[CH]C(F)(F)F(12166)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {10,S}
10 C u1 p0 c0 {9,S} {12,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-808.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,206,336,431,607,515,611,528,696,1312,1446,3000,3100,440,815,1455,1000,180,870.621],'cm^-1')),
        HinderedRotor(inertia=(0.185389,'amu*angstrom^2'), symmetry=1, barrier=(4.26245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00788716,'amu*angstrom^2'), symmetry=1, barrier=(34.5486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10721,'amu*angstrom^2'), symmetry=1, barrier=(48.449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.37287,0.0832398,-9.55809e-05,5.63995e-08,-1.33964e-11,-97137.8,27.6988], Tmin=(100,'K'), Tmax=(1016.25,'K')), NASAPolynomial(coeffs=[14.4521,0.0278242,-1.3788e-05,2.7437e-09,-1.97142e-13,-99999.5,-40.4527], Tmin=(1016.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-808.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Allyl_S) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = '[CH]C(F)=CC(F)C(F)(F)F(12197)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8  C u0 p0 c0 {6,S} {9,D} {12,S}
9  C u0 p0 c0 {5,S} {8,D} {10,S}
10 C u2 p0 c0 {9,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-719.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121836,0.0924304,-0.000120389,8.84655e-08,-2.68993e-11,-86458.3,29.0037], Tmin=(100,'K'), Tmax=(794.487,'K')), NASAPolynomial(coeffs=[10.7312,0.0390161,-1.95435e-05,3.84616e-09,-2.72707e-13,-88144.2,-19.7399], Tmin=(794.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCsCdF) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'F[CH]C(F)=C[CH]C(F)(F)F(12216)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,D} {12,S}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 C u1 p0 c0 {5,S} {9,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-798.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516089,0.0805238,-9.07743e-05,5.32749e-08,-1.26853e-11,-95868.2,27.126], Tmin=(100,'K'), Tmax=(1009.83,'K')), NASAPolynomial(coeffs=[13.4925,0.0291223,-1.44214e-05,2.86764e-09,-2.05899e-13,-98488.9,-35.6046], Tmin=(1009.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-798.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + radical(Allyl_S) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = '[CH]=C(F)C(F)C(F)[C](F)F(12217)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u1 p0 c0 {8,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-567.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,259,529,569,1128,1321,1390,3140,246,474,533,1155,190,488,555,1236,1407,3120,650,792.5,1650,180,180,1935.49],'cm^-1')),
        HinderedRotor(inertia=(0.403463,'amu*angstrom^2'), symmetry=1, barrier=(9.27641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402936,'amu*angstrom^2'), symmetry=1, barrier=(9.2643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81662,'amu*angstrom^2'), symmetry=1, barrier=(41.7677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.277788,0.103206,-0.000168814,1.46831e-07,-5.04764e-11,-68053.1,32.1886], Tmin=(100,'K'), Tmax=(799.658,'K')), NASAPolynomial(coeffs=[11.4092,0.0336739,-1.76175e-05,3.46532e-09,-2.42534e-13,-69568.2,-19.3681], Tmin=(799.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-567.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCsCdF) + group(Cds-CdsHH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'F[CH][C](F)C(F)C=C(F)F(3554)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u1 p0 c0 {2,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {10,D} {12,S}
9  C u1 p0 c0 {3,S} {7,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-622.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,212,367,445,1450,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,182,240,577,636,1210,1413,180,2046.01],'cm^-1')),
        HinderedRotor(inertia=(0.350555,'amu*angstrom^2'), symmetry=1, barrier=(8.05996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350386,'amu*angstrom^2'), symmetry=1, barrier=(8.05607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52917,'amu*angstrom^2'), symmetry=1, barrier=(35.1586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.253423,0.105359,-0.000182563,1.6729e-07,-5.96096e-11,-74703.4,33.2426], Tmin=(100,'K'), Tmax=(819.222,'K')), NASAPolynomial(coeffs=[9.37092,0.0374933,-2.00827e-05,3.96526e-09,-2.7712e-13,-75579.9,-6.99529], Tmin=(819.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = '[CH]=[C]C(F)C(F)C(F)(F)F(12218)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {10,D}
10 C u1 p0 c0 {9,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-584.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,164,312,561,654,898,1207,1299,3167,193,295,551,588,656,1146,1192,1350,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.383483,'amu*angstrom^2'), symmetry=1, barrier=(8.81703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383138,'amu*angstrom^2'), symmetry=1, barrier=(8.8091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.908083,'amu*angstrom^2'), symmetry=1, barrier=(20.8786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0363775,0.0988081,-0.000143702,1.00718e-07,-2.34126e-11,-70152.2,30.406], Tmin=(100,'K'), Tmax=(616.394,'K')), NASAPolynomial(coeffs=[12.7558,0.0311186,-1.62688e-05,3.22198e-09,-2.27429e-13,-72020.3,-27.4803], Tmin=(616.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-584.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sCsH)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = 'FC=[C]C(F)[CH]C(F)(F)F(12189)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-641.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,615,860,1140,1343,3152,1685,370,180,180,1978.46],'cm^-1')),
        HinderedRotor(inertia=(0.194081,'amu*angstrom^2'), symmetry=1, barrier=(4.4623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194243,'amu*angstrom^2'), symmetry=1, barrier=(4.46604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193904,'amu*angstrom^2'), symmetry=1, barrier=(4.45823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.227618,0.0930141,-0.00015318,1.4109e-07,-5.1694e-11,-77033.2,32.6662], Tmin=(100,'K'), Tmax=(787.305,'K')), NASAPolynomial(coeffs=[7.79389,0.0391004,-2.09831e-05,4.18878e-09,-2.96193e-13,-77745.1,1.0181], Tmin=(787.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-641.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cs_S) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH))"""),
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
    E0 = (-257.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-97.9873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (65.4097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-249.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-194.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (0.631511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-78.8675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-100.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (52.4187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-13.8728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-97.7769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (154.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (158.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-22.5692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-129.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-65.6617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-111.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-100.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-35.6971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-49.6399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-0.314684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-74.9862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['C2HF(58)', 'FC=CC(F)(F)F(1027)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(F)C([CH]F)C(F)(F)F(6905)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(F)(F)F(462)', '[CH]C(F)=CF(878)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['FC1=CC(C(F)(F)F)C1F(12175)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['C=C(F)C(F)=CC(F)(F)F(12209)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[CH]=C(F)C=CC(F)(F)F(12210)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]F(489)', 'FC=CC(F)(F)F(1027)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(161.758,'m^3/(mol*s)'), n=1.13767, Ea=(16.9891,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.011839139166081196, var=0.8671401158115553, Tref=1000.0, N=160, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', '[CH]=C(F)C(F)=CC(F)(F)F(12211)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.80541,'m^3/(mol*s)'), n=1.92214, Ea=(2.43569,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.14695202620752562, var=0.19149176908863236, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_Sp-7R!H=4CCClCl',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_Sp-7R!H=4CCClCl"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[CH]=C(F)C(F)C=C(F)F(9342)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.7375,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'C#CC(F)[CH]C(F)(F)F(12212)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(65.1617,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C2HF(58)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(2.70826,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]F(489)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -20.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F2(78)', '[CH]=C=C[CH]C(F)(F)F(12213)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(20.0613,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', '[CH]=C=C(F)[CH]C(F)(F)F(12214)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(295.678,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['[CH]C(F)=C(F)CC(F)(F)F(12215)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['[CH2]C(F)=C(F)[CH]C(F)(F)F(12166)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['[CH]C(F)=CC(F)C(F)(F)F(12197)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(146.543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['F[CH]C(F)=C[CH]C(F)(F)F(12216)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(157.23,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(F)C(F)C(F)[C](F)F(12217)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(146.025,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH][C](F)C(F)C=C(F)F(3554)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(187.347,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]C(F)C(F)C(F)(F)F(12218)'],
    products = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(198.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(F)C(F)[CH]C(F)(F)F(12108)'],
    products = ['FC=[C]C(F)[CH]C(F)(F)F(12189)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(182.804,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3642',
    isomers = [
        '[CH]=C(F)C(F)[CH]C(F)(F)F(12108)',
    ],
    reactants = [
        ('C2HF(58)', 'FC=CC(F)(F)F(1027)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3642',
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

