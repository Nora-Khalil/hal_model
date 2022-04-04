species(
    label = '[CH2]C(=C)C[C]=CF(8585)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,D} {14,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (247.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,615,860,1140,1343,3152,1685,370,215.419,980.336,980.338],'cm^-1')),
        HinderedRotor(inertia=(0.219486,'amu*angstrom^2'), symmetry=1, barrier=(7.22785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42326,'amu*angstrom^2'), symmetry=1, barrier=(13.9376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910886,'amu*angstrom^2'), symmetry=1, barrier=(29.996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93517,0.0615177,-4.97245e-05,2.08437e-08,-3.53692e-12,29851.8,25.5011], Tmin=(100,'K'), Tmax=(1395.86,'K')), NASAPolynomial(coeffs=[13.9083,0.0243416,-9.77452e-06,1.76331e-09,-1.19594e-13,26230.1,-41.4132], Tmin=(1395.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CF(5941)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (9.45959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1537.31],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2831,'J/mol'), sigma=(4.74838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=442.20 K, Pc=60 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96414,0.0021816,6.68618e-05,-1.27032e-07,7.57672e-11,1138.94,8.06307], Tmin=(10,'K'), Tmax=(518.59,'K')), NASAPolynomial(coeffs=[2.17835,0.0231528,-1.46137e-05,4.46872e-09,-5.27221e-13,1227.38,14.5728], Tmin=(518.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(9.45959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C=C(6556)',
    structure = adjacencyList("""1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37446,0.00704632,2.78302e-05,-3.99438e-08,1.55726e-11,21188.6,7.62051], Tmin=(100,'K'), Tmax=(949.709,'K')), NASAPolynomial(coeffs=[6.79959,0.00959973,-3.02065e-06,5.37819e-10,-3.92599e-14,19772.3,-12.7584], Tmin=(949.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(=C)C([CH2])=CF(8583)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {6,D}
3  C u0 p0 c0 {2,S} {5,S} {7,D}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {2,D} {8,S}
7  C u0 p0 c0 {3,D} {13,S} {14,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (129.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.591463,0.0631234,-3.33093e-05,-9.96746e-09,1.06108e-11,15762.5,22.1529], Tmin=(100,'K'), Tmax=(961.063,'K')), NASAPolynomial(coeffs=[18.1845,0.0176895,-5.76988e-06,1.01565e-09,-7.25536e-14,11097.5,-68.7012], Tmin=(961.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCFH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]CC[C]=CF(8588)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,D} {12,S}
5  C u0 p0 c0 {7,D} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {4,D}
7  C u1 p0 c0 {3,S} {5,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (356.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,615,860,1140,1343,3152,2950,3100,1380,975,1025,1650,1670,1700,300,440,293,294.343,295.092,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0592411,'amu*angstrom^2'), symmetry=1, barrier=(3.59713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43118,'amu*angstrom^2'), symmetry=1, barrier=(26.8623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304955,'amu*angstrom^2'), symmetry=1, barrier=(18.5213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3282.16,'J/mol'), sigma=(5.60824,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=512.67 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44414,0.0553973,-4.0543e-05,1.55334e-08,-2.49184e-12,42967.5,27.5036], Tmin=(100,'K'), Tmax=(1416.08,'K')), NASAPolynomial(coeffs=[10.8199,0.0289135,-1.24896e-05,2.32625e-09,-1.60182e-13,40312.1,-20.9906], Tmin=(1416.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-CsHH)(Cd-F1sH)) + radical(Cds_S)"""),
)

species(
    label = '[C]=CF(1334)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u2 p0 c0 {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (404.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([194,682,905,1196,1383,3221],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62308,0.0070472,-1.17409e-06,-1.98501e-09,8.12281e-13,48709.9,8.54797], Tmin=(100,'K'), Tmax=(1284.59,'K')), NASAPolynomial(coeffs=[5.40185,0.00468,-2.11337e-06,4.24439e-10,-3.06832e-14,47991.3,-1.49755], Tmin=(1284.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])=C(6712)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {4,D}
2  C u1 p0 c0 {1,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u0 p0 c0 {1,D} {5,S} {6,S}
5  H u0 p0 c0 {4,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (269.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.256303,'amu*angstrom^2'), symmetry=1, barrier=(104.754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0602849,'amu*angstrom^2'), symmetry=1, barrier=(24.6451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55518,0.0212333,2.47962e-05,-4.79596e-08,2.00952e-11,32418.8,14.0196], Tmin=(100,'K'), Tmax=(958.916,'K')), NASAPolynomial(coeffs=[10.7099,0.0134903,-4.19016e-06,7.65451e-10,-5.71873e-14,29646.9,-31.2785], Tmin=(958.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C]C[C]=CF(6834)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {6,D} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,D} {9,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {2,S} {3,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (372.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,615,860,1140,1343,3152,1670,1700,300,440,310.066,312.453,317.937],'cm^-1')),
        HinderedRotor(inertia=(0.0928139,'amu*angstrom^2'), symmetry=1, barrier=(6.32657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895553,'amu*angstrom^2'), symmetry=1, barrier=(6.33696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0622,0.0448898,-3.9659e-05,2.01535e-08,-4.42268e-12,44886.4,21.5872], Tmin=(100,'K'), Tmax=(1054.8,'K')), NASAPolynomial(coeffs=[7.3524,0.0248279,-1.11292e-05,2.12134e-09,-1.48779e-13,43770.4,-4.21715], Tmin=(1054.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C1CC(=CF)C1(9567)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u0 p0 c0 {1,S} {4,D} {12,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (4.34523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91662,0.0257397,6.31279e-05,-1.01685e-07,4.09081e-11,615.794,19.9693], Tmin=(100,'K'), Tmax=(960.322,'K')), NASAPolynomial(coeffs=[15.3805,0.0199934,-6.51787e-06,1.24347e-09,-9.59727e-14,-4291.11,-56.5255], Tmin=(960.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.34523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(CdCFH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C(C)C=C=CF(9380)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {7,D} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,D} {14,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (23.1521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.740292,0.0632717,-5.28364e-05,2.2739e-08,-3.90803e-12,2908.92,22.6723], Tmin=(100,'K'), Tmax=(1397.05,'K')), NASAPolynomial(coeffs=[15.2624,0.021693,-8.19432e-06,1.43628e-09,-9.59868e-14,-1148.74,-52.2441], Tmin=(1397.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.1521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'FC=[C]C[C]1CC1(9580)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {3,S} {4,S}
6  C u0 p0 c0 {1,S} {7,D} {14,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (329.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60568,0.0454901,-1.81148e-05,-2.32449e-09,2.39482e-12,39762.3,26.3987], Tmin=(100,'K'), Tmax=(1219.61,'K')), NASAPolynomial(coeffs=[10.2828,0.0287667,-1.198e-05,2.21157e-09,-1.52206e-13,36773,-20.7643], Tmin=(1219.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cyclopropane) + radical(Tertalkyl) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = 'C=C1C[C]([CH]F)C1(9577)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {3,S} {6,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {1,S} {4,S} {12,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (272.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24271,0.0444805,-2.08037e-05,3.6725e-09,-1.68172e-13,32853.5,20.6435], Tmin=(100,'K'), Tmax=(2199.32,'K')), NASAPolynomial(coeffs=[20.0621,0.0193921,-8.68558e-06,1.51267e-09,-9.46983e-14,23244.9,-83.3944], Tmin=(2199.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = '[CH2]C1([CH2])CC1=CF(9581)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u1 p0 c0 {2,S} {10,S} {11,S}
6  C u1 p0 c0 {2,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {4,D} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (351.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22799,0.0499576,-1.27225e-05,-1.90963e-08,1.06603e-11,42444.4,26.4788], Tmin=(100,'K'), Tmax=(1017.99,'K')), NASAPolynomial(coeffs=[13.9852,0.0237406,-9.32305e-06,1.74959e-09,-1.25097e-13,38608.1,-41.3797], Tmin=(1017.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cd(Cd-F)-Cs-Cs) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2][C]=C(6610)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {3,S} {6,S} {7,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (395.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.110742,'amu*angstrom^2'), symmetry=1, barrier=(23.6192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2737.61,'J/mol'), sigma=(4.88999,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=427.61 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28577,0.0103677,1.67968e-05,-2.70052e-08,1.05232e-11,47575.3,10.4038], Tmin=(100,'K'), Tmax=(980.89,'K')), NASAPolynomial(coeffs=[6.52797,0.0104288,-3.60851e-06,6.6848e-10,-4.85522e-14,46300.3,-8.4326], Tmin=(980.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,25474.2,-0.444973], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(=C)C=C=CF(9582)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {5,D}
3  C u0 p0 c0 {2,S} {7,D} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {2,D} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {7,D} {13,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (174.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,113,247,382,1207,3490,540,610,2055,515.639,515.754],'cm^-1')),
        HinderedRotor(inertia=(0.141127,'amu*angstrom^2'), symmetry=1, barrier=(26.6399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141262,'amu*angstrom^2'), symmetry=1, barrier=(26.6209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.893063,0.0591,-3.90423e-05,3.52419e-09,4.1787e-12,21125.7,21.7034], Tmin=(100,'K'), Tmax=(996.591,'K')), NASAPolynomial(coeffs=[15.7931,0.0184684,-6.74349e-06,1.2219e-09,-8.62366e-14,17203.8,-54.9067], Tmin=(996.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C][CH]F(5942)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (196.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.900551,'amu*angstrom^2'), symmetry=1, barrier=(20.7054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90543,0.0263147,-2.79302e-05,1.96926e-08,-6.19852e-12,23625.6,12.182], Tmin=(100,'K'), Tmax=(747.897,'K')), NASAPolynomial(coeffs=[4.81186,0.0161179,-7.47821e-06,1.46101e-09,-1.03887e-13,23340.4,3.53844], Tmin=(747.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'C#CCC([CH2])=C(7003)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {6,T}
6  C u0 p0 c0 {5,T} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (370.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,396.179],'cm^-1')),
        HinderedRotor(inertia=(0.204067,'amu*angstrom^2'), symmetry=1, barrier=(22.7334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.908337,'amu*angstrom^2'), symmetry=1, barrier=(101.158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204076,'amu*angstrom^2'), symmetry=1, barrier=(22.7333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41407,0.0486927,-2.27285e-05,-5.73777e-09,5.87842e-12,44670,19.6161], Tmin=(100,'K'), Tmax=(1006.52,'K')), NASAPolynomial(coeffs=[12.485,0.0217932,-8.12029e-06,1.46283e-09,-1.01801e-13,41575.4,-38.1686], Tmin=(1006.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)CC#CF(8633)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {10,S} {13,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {1,S} {6,T}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (270.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2175,525,239,401,1367,281.999,598.502],'cm^-1')),
        HinderedRotor(inertia=(7.30017e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868588,'amu*angstrom^2'), symmetry=1, barrier=(19.9711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.561492,'amu*angstrom^2'), symmetry=1, barrier=(31.6196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3894,0.0553427,-4.36076e-05,1.78916e-08,-3.02656e-12,32666.2,22.0267], Tmin=(100,'K'), Tmax=(1372.89,'K')), NASAPolynomial(coeffs=[11.6406,0.025476,-1.09764e-05,2.04651e-09,-1.41268e-13,29851.4,-30.6782], Tmin=(1372.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(CtCF) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])=CC=CF(9583)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  C u0 p0 c0 {3,D} {5,S} {6,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {7,D} {8,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {4,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (63.1612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32171,0.0462078,4.91752e-06,-4.22481e-08,2.06133e-11,7704.52,24.1411], Tmin=(100,'K'), Tmax=(947.25,'K')), NASAPolynomial(coeffs=[14.3563,0.0225208,-7.22457e-06,1.24153e-09,-8.70324e-14,3828.43,-45.4619], Tmin=(947.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.1612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(=C)CC=[C]F(9584)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {6,D}
4  C u0 p0 c0 {2,S} {7,D} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u0 p0 c0 {3,D} {13,S} {14,S}
7  C u1 p0 c0 {1,S} {4,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (259.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,167,640,1190,180.015,927.051],'cm^-1')),
        HinderedRotor(inertia=(0.267442,'amu*angstrom^2'), symmetry=1, barrier=(6.42328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.746215,'amu*angstrom^2'), symmetry=1, barrier=(17.9601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0548806,'amu*angstrom^2'), symmetry=1, barrier=(25.573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.855463,0.060878,-4.58421e-05,1.59626e-08,-1.74308e-12,31333.3,25.8051], Tmin=(100,'K'), Tmax=(1160.86,'K')), NASAPolynomial(coeffs=[14.3375,0.023489,-9.24534e-06,1.67351e-09,-1.14748e-13,27592.3,-43.8801], Tmin=(1160.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Allyl_P) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH2]C(C)=C[C]=CF(9384)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {7,S} {11,S}
5  C u1 p0 c0 {3,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,D} {14,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (168.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.942686,0.0611487,-4.98831e-05,2.20136e-08,-3.96236e-12,20421.8,25.1717], Tmin=(100,'K'), Tmax=(1322.13,'K')), NASAPolynomial(coeffs=[12.722,0.0255114,-9.45128e-06,1.62633e-09,-1.07346e-13,17307.1,-34.9457], Tmin=(1322.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(Cdj(Cd-CdH)(Cd-F1sH))"""),
)

species(
    label = '[CH]=C(C)C[C]=CF(8574)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {3,S} {7,D}
5  C u0 p0 c0 {1,S} {6,D} {13,S}
6  C u1 p0 c0 {2,S} {5,D}
7  C u1 p0 c0 {4,D} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (342.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,615,860,1140,1343,3152,1685,370,3120,650,792.5,1650,298.192,298.99],'cm^-1')),
        HinderedRotor(inertia=(0.135141,'amu*angstrom^2'), symmetry=1, barrier=(8.5371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134786,'amu*angstrom^2'), symmetry=1, barrier=(8.54223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134986,'amu*angstrom^2'), symmetry=1, barrier=(8.53693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3303.9,'J/mol'), sigma=(5.63277,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=516.06 K, Pc=41.95 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24378,0.0630148,-5.91842e-05,3.14447e-08,-7.05864e-12,41331.6,25.2124], Tmin=(100,'K'), Tmax=(1047.6,'K')), NASAPolynomial(coeffs=[9.52234,0.0314047,-1.39227e-05,2.64089e-09,-1.84782e-13,39597.1,-15.1116], Tmin=(1047.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)CC=CF(9385)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  C u2 p0 c0 {3,S} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (228.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848693,0.0606697,-3.74241e-05,9.51054e-09,-3.36123e-13,27613.4,26.8645], Tmin=(100,'K'), Tmax=(1223.41,'K')), NASAPolynomial(coeffs=[12.9574,0.0306363,-1.2318e-05,2.21466e-09,-1.50003e-13,23935.4,-36.9178], Tmin=(1223.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(C)C[C]=[C]F(9386)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {7,D}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (345.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,167,640,1190,340.838,340.846],'cm^-1')),
        HinderedRotor(inertia=(0.0964042,'amu*angstrom^2'), symmetry=1, barrier=(7.94763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964061,'amu*angstrom^2'), symmetry=1, barrier=(7.94769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964117,'amu*angstrom^2'), symmetry=1, barrier=(7.94766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31515,0.0621258,-5.86529e-05,3.2029e-08,-7.47639e-12,41691.9,25.0464], Tmin=(100,'K'), Tmax=(1004.37,'K')), NASAPolynomial(coeffs=[8.72505,0.0326149,-1.45787e-05,2.77383e-09,-1.94349e-13,40203.5,-10.7343], Tmin=(1004.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=[C]CC(=C)CF(9579)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u1 p0 c0 {2,S} {7,D}
7  C u1 p0 c0 {6,D} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (359.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,379.796,379.796],'cm^-1')),
        HinderedRotor(inertia=(0.115406,'amu*angstrom^2'), symmetry=1, barrier=(11.8129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115406,'amu*angstrom^2'), symmetry=1, barrier=(11.8129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115406,'amu*angstrom^2'), symmetry=1, barrier=(11.8129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34668,0.0616205,-5.54286e-05,2.795e-08,-6.02378e-12,43286.5,25.6379], Tmin=(100,'K'), Tmax=(1076.81,'K')), NASAPolynomial(coeffs=[9.17235,0.0325508,-1.49345e-05,2.87963e-09,-2.03284e-13,41601.1,-12.6956], Tmin=(1076.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)CC([CH2])=C(9585)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {6,D}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {12,S} {13,S}
7  C u1 p0 c0 {4,D} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (250.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,246,474,533,1155,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,188.366,1219.76],'cm^-1')),
        HinderedRotor(inertia=(0.00836951,'amu*angstrom^2'), symmetry=1, barrier=(8.82926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955106,'amu*angstrom^2'), symmetry=1, barrier=(23.7259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.969622,'amu*angstrom^2'), symmetry=1, barrier=(23.7298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671468,0.0642957,-5.3487e-05,2.26769e-08,-3.83041e-12,30310.3,25.8467], Tmin=(100,'K'), Tmax=(1421,'K')), NASAPolynomial(coeffs=[15.8879,0.0214621,-8.27153e-06,1.46363e-09,-9.8266e-14,25985.8,-52.91], Tmin=(1421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    E0 = (-3.12884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (166.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (276.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (423.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (503.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (5.15548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (75.1182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (228.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (123.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (101.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (161.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (147.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (153.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (252.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (232.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (341.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (153.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (178.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (199.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (221.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (391.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (128.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (255.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (180.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['C=C=CF(5941)', 'C=C=C(6556)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['[CH2]C(=C)C([CH2])=CF(8583)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]CC[C]=CF(8588)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=CF(1334)', '[CH2]C([CH2])=C(6712)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.13484e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(18)', 'C=[C]C[C]=CF(6834)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['C=C1CC(=CF)C1(9567)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['C=C(C)C=C=CF(9380)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.01596e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['FC=[C]C[C]1CC1(9580)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['C=C1C[C]([CH]F)C1(9577)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.06838e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['[CH2]C1([CH2])CC1=CF(9581)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(104.748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 103.3 to 104.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C=CF(5941)', '[CH2][C]=C(6610)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.0154,'m^3/(mol*s)'), n=2.41, Ea=(6.89117,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_N-Sp-5R!H-2C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_N-Sp-5R!H-2C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[CH2]C(=C)C=C=CF(9582)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1306.26,'m^3/(mol*s)'), n=1.34286, Ea=(11.7012,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2935929756704599, var=0.2435564102543282, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C][CH]F(5942)', 'C=C=C(6556)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.02833e-10,'m^3/(mol*s)'), n=4.64908, Ea=(31.7429,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24368554126510997, var=4.000885130961903, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing_3R->C_N-Sp-5R!H-4R!H_Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing_3R->C_N-Sp-5R!H-4R!H_Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'C#CCC([CH2])=C(7003)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.8932,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[CH2]C(=C)CC#CF(8633)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][CH]F(5942)', '[CH2][C]=C(6610)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['[CH2]C([CH2])=CC=CF(9583)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)CC=[C]F(9584)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['[CH2]C(C)=C[C]=CF(9384)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.19788e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(C)C[C]=CF(8574)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)C[C]=CF(8585)'],
    products = ['[CH]C(=C)CC=CF(9385)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.60477e+12,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C(C)C[C]=[C]F(9386)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]CC(=C)CF(9579)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(146.662,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(F)CC([CH2])=C(9585)'],
    products = ['[CH2]C(=C)C[C]=CF(8585)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(180.277,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3317',
    isomers = [
        '[CH2]C(=C)C[C]=CF(8585)',
    ],
    reactants = [
        ('C=C=CF(5941)', 'C=C=C(6556)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3317',
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

