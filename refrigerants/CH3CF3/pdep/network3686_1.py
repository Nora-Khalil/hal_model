species(
    label = 'F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {8,S} {9,S} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1164.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,253,522,585,610,849,1160,1215,1402,350,440,435,1725,3025,407.5,1350,352.5,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,895.882],'cm^-1')),
        HinderedRotor(inertia=(0.114581,'amu*angstrom^2'), symmetry=1, barrier=(2.63445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114735,'amu*angstrom^2'), symmetry=1, barrier=(2.63797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116454,'amu*angstrom^2'), symmetry=1, barrier=(2.67751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115786,'amu*angstrom^2'), symmetry=1, barrier=(2.66215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24489,0.126718,-0.0002025,1.77132e-07,-6.23565e-11,-139915,41.1412], Tmin=(100,'K'), Tmax=(758.85,'K')), NASAPolynomial(coeffs=[12.3064,0.0460503,-2.47863e-05,4.96573e-09,-3.5254e-13,-141705,-18.7442], Tmin=(758.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1164.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Cs_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[CH]C(=C(F)F)C([CH]F)C(F)(F)F(7322)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {4,S} {8,S} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1156.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,350,440,435,1725,334,575,1197,1424,3202,234,589,736,816,1240,3237,182,240,577,636,1210,1413,335.16,336.387,2210.62],'cm^-1')),
        HinderedRotor(inertia=(0.00328649,'amu*angstrom^2'), symmetry=1, barrier=(11.3929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142497,'amu*angstrom^2'), symmetry=1, barrier=(11.3852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486483,'amu*angstrom^2'), symmetry=1, barrier=(38.8741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.484686,'amu*angstrom^2'), symmetry=1, barrier=(38.8792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3417.16,'J/mol'), sigma=(5.74305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=533.75 K, Pc=40.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32446,0.126369,-0.000183957,1.40254e-07,-4.27733e-11,-138922,40.5802], Tmin=(100,'K'), Tmax=(801.733,'K')), NASAPolynomial(coeffs=[16.2291,0.0387894,-2.00967e-05,3.9965e-09,-2.84126e-13,-141736,-40.2264], Tmin=(801.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1156.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFF) + group(CsCsFHH) + group(CsCFHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(CdCFF) + radical(CsCsF1sH) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'FC(F)=[C]C(F)C(F)[CH]C(F)(F)F(12466)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {8,S} {13,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1046.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,180,858.696,858.9],'cm^-1')),
        HinderedRotor(inertia=(0.177888,'amu*angstrom^2'), symmetry=1, barrier=(4.08998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177894,'amu*angstrom^2'), symmetry=1, barrier=(4.09014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177692,'amu*angstrom^2'), symmetry=1, barrier=(4.0855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177695,'amu*angstrom^2'), symmetry=1, barrier=(4.08556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3418.75,'J/mol'), sigma=(5.75129,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=534.00 K, Pc=40.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.326448,0.110267,-0.000136675,6.55098e-08,7.45788e-12,-125719,37.8712], Tmin=(100,'K'), Tmax=(534.053,'K')), NASAPolynomial(coeffs=[10.7274,0.0493701,-2.71272e-05,5.52566e-09,-3.97867e-13,-127212,-11.4467], Tmin=(534.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1046.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH]C(=C[CH]C(F)(F)F)C(F)(F)F(12608)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {8,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {9,S} {12,S} {14,S}
12 C u0 p0 c0 {10,D} {11,S} {15,S}
13 C u1 p0 c0 {7,S} {10,S} {16,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1280.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.593598,0.108204,-0.000129674,8.12522e-08,-2.0685e-11,-153809,35.5533], Tmin=(100,'K'), Tmax=(946.615,'K')), NASAPolynomial(coeffs=[15.7818,0.0390077,-2.00265e-05,4.03123e-09,-2.90969e-13,-156910,-42.5502], Tmin=(946.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1280.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(CsCdFFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[CH]C([CH]F)=C(F)F(3373)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {8,D}
6  C u1 p0 c0 {1,S} {5,S} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-495.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,173,295,515,663,653,819,693,939,1188,1292,3205,3269,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.282253,'amu*angstrom^2'), symmetry=1, barrier=(6.48955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701374,'amu*angstrom^2'), symmetry=1, barrier=(56.9043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988945,0.0715839,-0.000105937,8.38657e-08,-2.66166e-11,-59462.3,23.9135], Tmin=(100,'K'), Tmax=(770.911,'K')), NASAPolynomial(coeffs=[10.3313,0.0231078,-1.16114e-05,2.29254e-09,-1.62266e-13,-60902.7,-18.7271], Tmin=(770.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-495.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'FC(F)=[C]C(F)[CH]C(F)(F)F(3367)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {6,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-832.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,1480.67,1480.96],'cm^-1')),
        HinderedRotor(inertia=(0.271278,'amu*angstrom^2'), symmetry=1, barrier=(6.23721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271037,'amu*angstrom^2'), symmetry=1, barrier=(6.23167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271168,'amu*angstrom^2'), symmetry=1, barrier=(6.23469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0915714,0.0960851,-0.000157193,1.45e-07,-5.38228e-11,-100048,35.2443], Tmin=(100,'K'), Tmax=(761.998,'K')), NASAPolynomial(coeffs=[8.12054,0.0408052,-2.25211e-05,4.55864e-09,-3.25455e-13,-100890,1.19334], Tmin=(761.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-832.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFF) + radical(Cs_S) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC(F)=C1C(F)C(C(F)(F)F)C1F(12556)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {12,S} {15,S}
10 C u0 p0 c0 {2,S} {8,S} {12,S} {16,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
12 C u0 p0 c0 {9,S} {10,S} {13,D}
13 C u0 p0 c0 {6,S} {7,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1367.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.174427,0.102098,-0.000120621,7.88683e-08,-2.17384e-11,-164293,28.7369], Tmin=(100,'K'), Tmax=(862.482,'K')), NASAPolynomial(coeffs=[11.6454,0.0472787,-2.52796e-05,5.17183e-09,-3.76295e-13,-166332,-26.5382], Tmin=(862.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1367.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCCFH) + group(CsCsFFF) + group(Cds-CdsCsCs) + group(CdCFF) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FCC(C(F)=CC(F)(F)F)=C(F)F(12609)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {4,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {11,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u0 p0 c0 {9,S} {12,D} {16,S}
12 C u0 p0 c0 {5,S} {10,S} {11,D}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1379.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.745822,0.11193,-0.000138944,8.95469e-08,-2.33214e-11,-165763,34.7824], Tmin=(100,'K'), Tmax=(928.607,'K')), NASAPolynomial(coeffs=[16.3803,0.038158,-1.97777e-05,3.99418e-09,-2.88662e-13,-168943,-46.5726], Tmin=(928.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1379.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF)"""),
)

species(
    label = 'FC([CH]C(F)(F)F)[C]1C(F)C1(F)F(12610)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {9,S} {12,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {12,S} {13,S} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {13,S}
12 C u1 p0 c0 {8,S} {9,S} {10,S}
13 C u1 p0 c0 {10,S} {11,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1077.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258285,0.100332,-8.38698e-05,-6.42249e-08,1.17586e-10,-129485,34.3928], Tmin=(100,'K'), Tmax=(469.122,'K')), NASAPolynomial(coeffs=[10.2389,0.048684,-2.5688e-05,5.09563e-09,-3.5918e-13,-130790,-10.1269], Tmin=(469.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1077.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH][C]1C(F)C(C(F)(F)F)C1(F)F(12611)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {12,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
12 C u1 p0 c0 {9,S} {10,S} {13,S}
13 C u1 p0 c0 {7,S} {12,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1159.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.53336,0.109429,-0.000148125,1.13268e-07,-3.59363e-11,-139328,35.1116], Tmin=(100,'K'), Tmax=(761.272,'K')), NASAPolynomial(coeffs=[11.423,0.0466123,-2.43635e-05,4.89746e-09,-3.51404e-13,-141149,-19.3111], Tmin=(761.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1159.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCsFHH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1([C](F)F)C(F)C1C(F)(F)F(12597)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {8,S} {10,S} {11,S} {14,S}
10 C u0 p0 c0 {1,S} {8,S} {9,S} {15,S}
11 C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
12 C u1 p0 c0 {5,S} {8,S} {16,S}
13 C u1 p0 c0 {6,S} {7,S} {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1109.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06098,0.119305,-0.000150537,9.65823e-08,-2.48739e-11,-133302,35.4869], Tmin=(100,'K'), Tmax=(940.909,'K')), NASAPolynomial(coeffs=[18.1997,0.037422,-1.99955e-05,4.08724e-09,-2.9733e-13,-136927,-56.2613], Tmin=(940.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1109.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFFF) + ring(Cs-Cs(C-FFF)-Cs) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[CH]C(C=CC(F)(F)F)=C(F)F(12612)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
8  C u0 p0 c0 {10,S} {11,S} {12,D}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 C u1 p0 c0 {4,S} {8,S} {15,S}
12 C u0 p0 c0 {5,S} {6,S} {8,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1050.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,234,589,736,816,1240,3237,182,240,577,636,1210,1413,204.973,205.024],'cm^-1')),
        HinderedRotor(inertia=(0.00402916,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22869,'amu*angstrom^2'), symmetry=1, barrier=(36.613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533602,'amu*angstrom^2'), symmetry=1, barrier=(15.8747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.647,0.106203,-0.000130156,8.06164e-08,-1.98446e-11,-126146,32.7792], Tmin=(100,'K'), Tmax=(988.613,'K')), NASAPolynomial(coeffs=[18.1023,0.0303404,-1.50497e-05,2.99327e-09,-2.14963e-13,-129854,-57.4601], Tmin=(988.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1050.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCdFFF) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + radical(CsCdF1sH)"""),
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
    label = 'F[CH]C(C(F)=CC(F)(F)F)=C(F)F(12613)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
9  C u0 p0 c0 {11,S} {12,S} {13,D}
10 C u0 p0 c0 {8,S} {11,D} {14,S}
11 C u0 p0 c0 {4,S} {9,S} {10,D}
12 C u1 p0 c0 {5,S} {9,S} {15,S}
13 C u0 p0 c0 {6,S} {7,S} {9,D}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1239.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,350,440,435,1725,3010,987.5,1337.5,450,1655,280,518,736,852,873,234,589,736,816,1240,3237,182,240,577,636,1210,1413,268.866,2573.66],'cm^-1')),
        HinderedRotor(inertia=(0.16064,'amu*angstrom^2'), symmetry=1, barrier=(7.97855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754244,'amu*angstrom^2'), symmetry=1, barrier=(35.1158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711153,'amu*angstrom^2'), symmetry=1, barrier=(35.0981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (207.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.80762,0.112378,-0.00014223,9.11965e-08,-2.33701e-11,-148866,34.8045], Tmin=(100,'K'), Tmax=(947.718,'K')), NASAPolynomial(coeffs=[17.8354,0.0336927,-1.769e-05,3.5895e-09,-2.60126e-13,-152399,-54.136], Tmin=(947.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1239.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF) + radical(CsCdF1sH)"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)C=C(F)F(9367)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {7,S} {12,D} {14,S}
10 C u1 p0 c0 {4,S} {8,S} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {8,D}
12 C u0 p0 c0 {2,S} {3,S} {9,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-968.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,141,223,164,316,539,615,578,694,1133,1287,1372,1454,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.907634,'amu*angstrom^2'), symmetry=1, barrier=(20.8683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17357,'amu*angstrom^2'), symmetry=1, barrier=(49.9746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0183755,'amu*angstrom^2'), symmetry=1, barrier=(5.5238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.385492,0.1081,-0.00014861,1.01934e-07,-2.30343e-11,-116302,34.8919], Tmin=(100,'K'), Tmax=(600.007,'K')), NASAPolynomial(coeffs=[11.8401,0.0420653,-2.21972e-05,4.44442e-09,-3.16889e-13,-118048,-20.1648], Tmin=(600.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-968.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[CH]C([C](F)CC(F)(F)F)=C(F)F(12614)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {11,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u0 p0 c0 {11,S} {12,S} {13,D}
11 C u1 p0 c0 {4,S} {8,S} {10,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1221.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10536,0.122046,-0.000170732,1.26594e-07,-3.79292e-11,-146724,38.2924], Tmin=(100,'K'), Tmax=(812.115,'K')), NASAPolynomial(coeffs=[15.2581,0.0414524,-2.18785e-05,4.40445e-09,-3.16066e-13,-149382,-37.247], Tmin=(812.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1221.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CsCdCsF1s) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'FCC([C](F)F)=C(F)[CH]C(F)(F)F(12615)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {9,S} {12,S} {16,S}
12 C u0 p0 c0 {5,S} {10,D} {11,S}
13 C u1 p0 c0 {6,S} {7,S} {10,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1207.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.574853,0.10966,-0.000136725,9.09831e-08,-2.48214e-11,-145120,37.7175], Tmin=(100,'K'), Tmax=(882.471,'K')), NASAPolynomial(coeffs=[14.2856,0.0423073,-2.22517e-05,4.51091e-09,-3.26344e-13,-147743,-32.1187], Tmin=(882.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1207.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C([CH]C(F)C(F)(F)F)=C(F)F(12591)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
10 C u0 p0 c0 {11,S} {12,S} {13,D}
11 C u1 p0 c0 {8,S} {10,S} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1202.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.922511,0.116196,-0.000150438,1.01751e-07,-2.77452e-11,-144463,36.5184], Tmin=(100,'K'), Tmax=(890.03,'K')), NASAPolynomial(coeffs=[16.3098,0.0387524,-1.9923e-05,3.99341e-09,-2.86933e-13,-147531,-44.6105], Tmin=(890.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1202.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C](F)C(=C[CH]C(F)(F)F)C(F)F(12616)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {9,S} {12,S} {15,S}
12 C u0 p0 c0 {10,D} {11,S} {16,S}
13 C u1 p0 c0 {6,S} {7,S} {10,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1224.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.54604,0.10961,-0.000138974,9.52876e-08,-2.69191e-11,-147061,37.7613], Tmin=(100,'K'), Tmax=(851.945,'K')), NASAPolynomial(coeffs=[13.4546,0.0438748,-2.32341e-05,4.71726e-09,-3.41319e-13,-149447,-27.54], Tmin=(851.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1224.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(CsCFFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)C(F)[C](F)F(12617)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {2,S} {9,S} {10,S} {15,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {3,S} {4,S} {9,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1088.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,259,529,569,1128,1321,1390,3140,350,440,435,1725,190,488,555,1236,1407,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,180,2376.73],'cm^-1')),
        HinderedRotor(inertia=(0.194866,'amu*angstrom^2'), symmetry=1, barrier=(4.48036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196298,'amu*angstrom^2'), symmetry=1, barrier=(4.51327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71024,'amu*angstrom^2'), symmetry=1, barrier=(39.3218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71192,'amu*angstrom^2'), symmetry=1, barrier=(39.3604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59595,0.135076,-0.000221169,1.93026e-07,-6.69345e-11,-130754,41.2717], Tmin=(100,'K'), Tmax=(782.311,'K')), NASAPolynomial(coeffs=[13.7137,0.0444796,-2.38437e-05,4.74381e-09,-3.34541e-13,-132772,-26.4209], Tmin=(782.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1088.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C](F)[C](C(F)F)C(F)C=C(F)F(7281)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {10,S} {15,S}
10 C u1 p0 c0 {8,S} {9,S} {12,S}
11 C u0 p0 c0 {8,S} {13,D} {16,S}
12 C u1 p0 c0 {4,S} {5,S} {10,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1102.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,522,611,926,1093,1137,1374,1416,3112,360,370,350,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,182,240,577,636,1210,1413,180,180,1778.76],'cm^-1')),
        HinderedRotor(inertia=(0.324074,'amu*angstrom^2'), symmetry=1, barrier=(7.45111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323826,'amu*angstrom^2'), symmetry=1, barrier=(7.4454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32398,'amu*angstrom^2'), symmetry=1, barrier=(7.44895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323944,'amu*angstrom^2'), symmetry=1, barrier=(7.4481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59478,0.141661,-0.000258398,2.43388e-07,-8.75837e-11,-132398,41.9531], Tmin=(100,'K'), Tmax=(837.048,'K')), NASAPolynomial(coeffs=[9.55383,0.0505575,-2.73518e-05,5.38258e-09,-3.7387e-13,-132939,-1.93353], Tmin=(837.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1102.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]C(=CF)C(F)C(F)C(F)(F)F(12618)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u0 p0 c0 {9,S} {12,D} {13,S}
12 C u0 p0 c0 {6,S} {11,D} {16,S}
13 C u2 p0 c0 {7,S} {11,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1103.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,284,328,853,1146,1135,1297,3239,193,295,551,588,656,1146,1192,1350,350,440,435,1725,194,682,905,1196,1383,3221,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66912,0.13638,-0.000223418,1.93177e-07,-6.63501e-11,-132502,37.7876], Tmin=(100,'K'), Tmax=(776.184,'K')), NASAPolynomial(coeffs=[14.6322,0.0424985,-2.29075e-05,4.56881e-09,-3.22684e-13,-134735,-34.8106], Tmin=(776.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1103.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=C(C(F)F)C(F)[CH]C(F)(F)F(12619)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {11,S} {12,S} {14,S}
9  C u0 p0 c0 {5,S} {6,S} {11,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {12,S}
11 C u0 p0 c0 {8,S} {9,S} {13,D}
12 C u1 p0 c0 {8,S} {10,S} {16,S}
13 C u1 p0 c0 {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1077.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,171,205,598,1104,1143,1317,1411,3153,253,522,585,610,849,1160,1215,1402,350,440,435,1725,3025,407.5,1350,352.5,167,640,1190,180,180,1603.18],'cm^-1')),
        HinderedRotor(inertia=(0.214466,'amu*angstrom^2'), symmetry=1, barrier=(4.93099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21049,'amu*angstrom^2'), symmetry=1, barrier=(4.83957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213871,'amu*angstrom^2'), symmetry=1, barrier=(4.91731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211423,'amu*angstrom^2'), symmetry=1, barrier=(4.86102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51398,0.136067,-0.000235481,2.15279e-07,-7.68953e-11,-129415,41.5482], Tmin=(100,'K'), Tmax=(807.508,'K')), NASAPolynomial(coeffs=[11.5019,0.0474133,-2.58845e-05,5.15839e-09,-3.62775e-13,-130729,-13.5817], Tmin=(807.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1077.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Cs_S) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)[CH]C(F)(F)F(12467)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {15,S}
12 C u0 p0 c0 {7,S} {13,D} {16,S}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1071.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,136,307,446,511,682,757,1180,1185,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,615,860,1140,1343,3152,1685,370,180,180,180,1992.07],'cm^-1')),
        HinderedRotor(inertia=(0.232135,'amu*angstrom^2'), symmetry=1, barrier=(5.33723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231431,'amu*angstrom^2'), symmetry=1, barrier=(5.32106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232178,'amu*angstrom^2'), symmetry=1, barrier=(5.33822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30284,'amu*angstrom^2'), symmetry=1, barrier=(29.9549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3407.02,'J/mol'), sigma=(6.03455,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=532.17 K, Pc=35.18 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38759,0.129761,-0.000207158,1.78751e-07,-6.21135e-11,-128689,39.8177], Tmin=(100,'K'), Tmax=(743.947,'K')), NASAPolynomial(coeffs=[13.4468,0.0445083,-2.41931e-05,4.86901e-09,-3.46712e-13,-130744,-26.3404], Tmin=(743.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1071.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
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
    label = 'FC=C1C(F)C(C(F)(F)F)C1(F)F(12524)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {8,S} {12,S} {15,S}
11 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
12 C u0 p0 c0 {9,S} {10,S} {13,D}
13 C u0 p0 c0 {7,S} {12,D} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1381.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.279315,0.103362,-0.000118953,7.31239e-08,-1.86633e-11,-165998,29.5383], Tmin=(100,'K'), Tmax=(933.187,'K')), NASAPolynomial(coeffs=[13.6124,0.0438173,-2.32412e-05,4.7477e-09,-3.45442e-13,-168590,-36.5207], Tmin=(933.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1381.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFF) + group(CsCCFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FC=C(C(F)=CC(F)(F)F)C(F)F(12620)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {4,S} {5,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {11,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u0 p0 c0 {9,S} {12,D} {15,S}
12 C u0 p0 c0 {6,S} {10,S} {11,D}
13 C u0 p0 c0 {7,S} {10,D} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1397.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.768288,0.112182,-0.000138852,8.88542e-08,-2.29372e-11,-167858,34.9697], Tmin=(100,'K'), Tmax=(937.221,'K')), NASAPolynomial(coeffs=[16.7085,0.0375905,-1.94689e-05,3.9324e-09,-2.84288e-13,-171134,-48.2122], Tmin=(937.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1397.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFH)"""),
)

species(
    label = 'F[C](F)[C]1C(F)C(C(F)(F)F)C1F(12621)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {12,S} {15,S}
10 C u0 p0 c0 {2,S} {8,S} {12,S} {16,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
12 C u1 p0 c0 {9,S} {10,S} {13,S}
13 C u1 p0 c0 {6,S} {7,S} {12,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1166.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.514567,0.108924,-0.000151106,1.1977e-07,-3.94206e-11,-140167,35.7525], Tmin=(100,'K'), Tmax=(734.971,'K')), NASAPolynomial(coeffs=[10.8863,0.046877,-2.44761e-05,4.90991e-09,-3.51503e-13,-141843,-15.7394], Tmin=(734.971,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1166.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsFFF) + group(CsCsFFH) + ring(Cs(F)-Cs-Cs-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
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
    label = 'F[CH]C(=C(F)[CH]C(F)(F)F)C(F)F(12622)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {14,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {9,S} {12,S} {15,S}
12 C u0 p0 c0 {6,S} {10,D} {11,S}
13 C u1 p0 c0 {7,S} {10,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1233.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.633655,0.110301,-0.000136941,8.98463e-08,-2.40324e-11,-148234,37.0641], Tmin=(100,'K'), Tmax=(901.637,'K')), NASAPolynomial(coeffs=[15.0446,0.0407465,-2.12283e-05,4.2892e-09,-3.09855e-13,-151061,-36.9512], Tmin=(901.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1233.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C]C(=C(F)F)C(F)CC(F)(F)F(12623)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {16,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
11 C u0 p0 c0 {9,S} {12,D} {13,S}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 C u2 p0 c0 {7,S} {11,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
"""),
    E0 = (-1134.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,284,328,853,1146,1135,1297,3239,193,295,551,588,656,1146,1192,1350,350,440,435,1725,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53016,0.134346,-0.000223953,1.98707e-07,-6.97646e-11,-136266,38.0716], Tmin=(100,'K'), Tmax=(788.614,'K')), NASAPolynomial(coeffs=[12.9862,0.0450763,-2.44079e-05,4.86992e-09,-3.43661e-13,-138069,-25.4305], Tmin=(788.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1134.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[CH][C](C(F)C=C(F)F)C(F)(F)F(12566)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {12,S}
11 C u0 p0 c0 {8,S} {13,D} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1133.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,253,522,585,610,849,1160,1215,1402,360,370,350,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,182,240,577,636,1210,1413,180,180,1958.72],'cm^-1')),
        HinderedRotor(inertia=(0.20524,'amu*angstrom^2'), symmetry=1, barrier=(4.71888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205227,'amu*angstrom^2'), symmetry=1, barrier=(4.71857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205227,'amu*angstrom^2'), symmetry=1, barrier=(4.71858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20522,'amu*angstrom^2'), symmetry=1, barrier=(4.71842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47624,0.139156,-0.000253687,2.40575e-07,-8.71736e-11,-136196,40.6148], Tmin=(100,'K'), Tmax=(835.539,'K')), NASAPolynomial(coeffs=[8.73412,0.0518702,-2.8041e-05,5.52236e-09,-3.8396e-13,-136562,1.21255], Tmin=(835.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1133.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = '[CH]C(=C(F)F)C(F)C(F)C(F)(F)F(12624)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u0 p0 c0 {9,S} {12,D} {13,S}
12 C u0 p0 c0 {6,S} {7,S} {11,D}
13 C u2 p0 c0 {11,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1128.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,284,328,853,1146,1135,1297,3239,193,295,551,588,656,1146,1192,1350,350,440,435,1725,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46922,0.130114,-0.000184655,1.42616e-07,-4.46937e-11,-135564,37.941], Tmin=(100,'K'), Tmax=(777.197,'K')), NASAPolynomial(coeffs=[14.7154,0.0468164,-2.38866e-05,4.71047e-09,-3.33287e-13,-138080,-36.0606], Tmin=(777.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1128.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C(F)[CH]C(F)(F)F)C(F)(F)F(12625)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  F u0 p3 c0 {9,S}
8  C u0 p0 c0 {1,S} {11,S} {12,S} {14,S}
9  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {12,S}
11 C u0 p0 c0 {8,S} {9,S} {13,D}
12 C u1 p0 c0 {8,S} {10,S} {15,S}
13 C u1 p0 c0 {11,D} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1137.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,219,296,586,564,718,793,1177,1228,253,522,585,610,849,1160,1215,1402,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,180,1462.32],'cm^-1')),
        HinderedRotor(inertia=(0.259755,'amu*angstrom^2'), symmetry=1, barrier=(5.97227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259466,'amu*angstrom^2'), symmetry=1, barrier=(5.96562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259347,'amu*angstrom^2'), symmetry=1, barrier=(5.9629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259536,'amu*angstrom^2'), symmetry=1, barrier=(5.96723,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4319,0.131131,-0.000214053,1.86948e-07,-6.50663e-11,-136661,39.9218], Tmin=(100,'K'), Tmax=(773.854,'K')), NASAPolynomial(coeffs=[13.3787,0.0437161,-2.3561e-05,4.70699e-09,-3.32991e-13,-138628,-25.6325], Tmin=(773.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1137.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFF) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cds_P)"""),
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
    E0 = (-347.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-179.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-134.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-94.0928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (2.67489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (199.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-339.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-284.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-116.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-222.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-282.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-97.5528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-209.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-205.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-37.4034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-128.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (35.3896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (123.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-219.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-138.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-183.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-149.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-125.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-117.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-74.3454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-74.4012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-128.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (209.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-339.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-284.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-222.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-28.0687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-147.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-273.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-134.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-101.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-106.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FC=C=C(F)F(1325)', 'FC=CC(F)(F)F(1027)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(=C(F)F)C([CH]F)C(F)(F)F(7322)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC(F)=[C]C(F)C(F)[CH]C(F)(F)F(12466)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[CH]C(=C[CH]C(F)(F)F)C(F)(F)F(12608)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72966e+11,'s^-1'), n=0.63878, Ea=(253.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(F)(F)F(462)', 'F[CH]C([CH]F)=C(F)F(3373)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(804)', 'FC(F)=[C]C(F)[CH]C(F)(F)F(3367)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FC(F)=C1C(F)C(C(F)(F)F)C1F(12556)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FCC(C(F)=CC(F)(F)F)=C(F)F(12609)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FC([CH]C(F)(F)F)[C]1C(F)C1(F)F(12610)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[CH][C]1C(F)C(C(F)(F)F)C1(F)F(12611)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[CH]C1([C](F)F)C(F)C1C(F)(F)F(12597)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(64.7047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[CH]C(C=CC(F)(F)F)=C(F)F(12612)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.5786,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C]=C(F)F(2842)', 'FC=CC(F)(F)F(1027)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(5.55744,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'F[CH]C(C(F)=CC(F)(F)F)=C(F)F(12613)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.80541,'m^3/(mol*s)'), n=1.92214, Ea=(4.25282,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.14695202620752562, var=0.19149176908863236, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_Sp-7R!H=4CCClCl',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_Sp-7R!H=4CCClCl"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'F[CH]C(=C(F)F)C(F)C=C(F)F(9367)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.7375,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['FC=C=C(F)F(1325)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -12.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHF(40)', 'FC(F)=[C]C(F)[CH]C(F)(F)F(3367)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[CH]C([C](F)CC(F)(F)F)=C(F)F(12614)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FCC([C](F)F)=C(F)[CH]C(F)(F)F(12615)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.18083e+09,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[CH]C([CH]C(F)C(F)(F)F)=C(F)F(12591)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(163.736,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[C](F)C(=C[CH]C(F)(F)F)C(F)F(12616)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(198.514,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(=C(F)F)C(F)C(F)[C](F)F(12617)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(146.025,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C](F)[C](C(F)F)C(F)C=C(F)F(7281)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(167.721,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C]C(=CF)C(F)C(F)C(F)(F)F(12618)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(211.759,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]=C(C(F)F)C(F)[CH]C(F)(F)F(12619)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C(F)[CH]C(F)(F)F(12467)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[C]F(138)', 'FC=[C]C(F)[CH]C(F)(F)F(12189)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FC=C1C(F)C(C(F)(F)F)C1(F)F(12524)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['FC=C(C(F)=CC(F)(F)F)C(F)F(12620)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[C](F)[C]1C(F)C(C(F)(F)F)C1F(12621)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CF2(43)', 'FC=[C]C(F)[CH]C(F)(F)F(12189)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    products = ['F[CH]C(=C(F)[CH]C(F)(F)F)C(F)F(12622)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.97418e+09,'s^-1'), n=1.23333, Ea=(200.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['F[C]C(=C(F)F)C(F)CC(F)(F)F(12623)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F[CH][C](C(F)C=C(F)F)C(F)(F)F(12566)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.94559e+07,'s^-1'), n=1.15307, Ea=(182.54,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]C(=C(F)F)C(F)C(F)C(F)(F)F(12624)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(209.724,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C(C(F)[CH]C(F)(F)F)C(F)(F)F(12625)'],
    products = ['F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #3686',
    isomers = [
        'F[CH]C(=C(F)F)C(F)[CH]C(F)(F)F(12465)',
    ],
    reactants = [
        ('FC=C=C(F)F(1325)', 'FC=CC(F)(F)F(1027)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3686',
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

