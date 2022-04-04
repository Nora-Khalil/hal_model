species(
    label = 'F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {13,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
11 C u1 p0 c0 {5,S} {8,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1056.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,193,295,551,588,656,1146,1192,1350,334,575,1197,1424,3202,562,600,623,1070,1265,1685,370,180,180,1633.04,1633.77],'cm^-1')),
        HinderedRotor(inertia=(0.214463,'amu*angstrom^2'), symmetry=1, barrier=(4.93092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214125,'amu*angstrom^2'), symmetry=1, barrier=(4.92315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214159,'amu*angstrom^2'), symmetry=1, barrier=(4.92394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64351,'amu*angstrom^2'), symmetry=1, barrier=(37.7875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18841,0.125862,-0.000202175,1.79179e-07,-6.39879e-11,-126860,40.9109], Tmin=(100,'K'), Tmax=(757.848,'K')), NASAPolynomial(coeffs=[11.5829,0.0476154,-2.58513e-05,5.19849e-09,-3.69869e-13,-128484,-15.1087], Tmin=(757.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1056.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
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
    label = 'F[CH]C(F)C(F)[C]=C(F)F(3030)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {12,S}
8  C u1 p0 c0 {3,S} {6,S} {13,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-549.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,334,575,1197,1424,3202,562,600,623,1070,1265,1685,370,180,180,1510.28],'cm^-1')),
        HinderedRotor(inertia=(0.268661,'amu*angstrom^2'), symmetry=1, barrier=(6.17704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26853,'amu*angstrom^2'), symmetry=1, barrier=(6.17403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67366,'amu*angstrom^2'), symmetry=1, barrier=(38.4807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3185.22,'J/mol'), sigma=(5.13817,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=497.52 K, Pc=53.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.652361,0.0862597,-0.000101708,3.22885e-08,2.30117e-11,-66005.2,30.5731], Tmin=(100,'K'), Tmax=(517.923,'K')), NASAPolynomial(coeffs=[9.59627,0.0376914,-2.04364e-05,4.12443e-09,-2.94761e-13,-67206.7,-9.34688], Tmin=(517.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-549.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.326448,0.110267,-0.000136675,6.55098e-08,7.45788e-12,-125719,37.8712], Tmin=(100,'K'), Tmax=(534.053,'K')), NASAPolynomial(coeffs=[10.7274,0.0493701,-2.71272e-05,5.52566e-09,-3.97867e-13,-127212,-11.4467], Tmin=(534.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1046.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC(F)=[C]C(F)[CH]C(F)C(F)(F)F(12555)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {11,S} {13,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u1 p0 c0 {8,S} {9,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1041.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,174,267,591,721,1107,1278,1348,3273,193,295,551,588,656,1146,1192,1350,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,180,1179.63,1179.64],'cm^-1')),
        HinderedRotor(inertia=(0.216784,'amu*angstrom^2'), symmetry=1, barrier=(4.98429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216735,'amu*angstrom^2'), symmetry=1, barrier=(4.98315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21692,'amu*angstrom^2'), symmetry=1, barrier=(4.98741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216691,'amu*angstrom^2'), symmetry=1, barrier=(4.98215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0214281,0.107271,-0.000104512,-3.63803e-08,1.00819e-10,-125099,38.0218], Tmin=(100,'K'), Tmax=(476.417,'K')), NASAPolynomial(coeffs=[11.1146,0.0480902,-2.62274e-05,5.2661e-09,-3.73532e-13,-126549,-11.5332], Tmin=(476.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1041.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCCFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(Cs_S) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32446,0.126369,-0.000183957,1.40254e-07,-4.27733e-11,-138922,40.5802], Tmin=(100,'K'), Tmax=(801.733,'K')), NASAPolynomial(coeffs=[16.2291,0.0387894,-2.00967e-05,3.9965e-09,-2.84126e-13,-141736,-40.2264], Tmin=(801.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1156.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFF) + group(CsCsFHH) + group(CsCFHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(CdCFF) + radical(CsCsF1sH) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = '[C]=C(F)F(1205)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u2 p0 c0 {3,D}
"""),
    E0 = (195.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([182,240,577,636,1210,1413],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28657,0.0148134,-1.42116e-05,6.41427e-09,-1.12277e-12,23593,10.1566], Tmin=(100,'K'), Tmax=(1382.13,'K')), NASAPolynomial(coeffs=[7.28631,0.0032378,-1.64877e-06,3.54597e-10,-2.66861e-14,22487.4,-10.4342], Tmin=(1382.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'F[CH]C([CH]F)C(F)(F)F(12266)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u1 p0 c0 {4,S} {6,S} {11,S}
9  C u1 p0 c0 {5,S} {6,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-790.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,262,406,528,622,1148,1246,1368,1480,3164,3240,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.245899,'amu*angstrom^2'), symmetry=1, barrier=(5.65371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88562,'amu*angstrom^2'), symmetry=1, barrier=(43.3541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00220658,'amu*angstrom^2'), symmetry=1, barrier=(5.65653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.112557,0.0980113,-0.000176273,1.66121e-07,-6.0163e-11,-95001.7,26.7801], Tmin=(100,'K'), Tmax=(826.425,'K')), NASAPolynomial(coeffs=[7.895,0.0359539,-1.93678e-05,3.83787e-09,-2.68378e-13,-95455.1,-4.24305], Tmin=(826.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-790.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    label = 'FC=C(C(F)C=C(F)F)C(F)(F)F(12557)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u0 p0 c0 {8,S} {13,D} {15,S}
12 C u0 p0 c0 {5,S} {10,D} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1380.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.864873,0.119669,-0.000168113,1.16355e-07,-2.65885e-11,-165825,35.1775], Tmin=(100,'K'), Tmax=(605.253,'K')), NASAPolynomial(coeffs=[13.3483,0.0433963,-2.28524e-05,4.56203e-09,-3.24332e-13,-167869,-28.9284], Tmin=(605.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1380.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF)"""),
)

species(
    label = 'FCC(C(F)=C=C(F)F)C(F)(F)F(12558)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
11 C u0 p0 c0 {5,S} {8,S} {13,D}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u0 p0 c0 {11,D} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1307.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.825122,0.112667,-0.000142352,9.31031e-08,-2.44363e-11,-157063,37.433], Tmin=(100,'K'), Tmax=(925.355,'K')), NASAPolynomial(coeffs=[16.8988,0.0360496,-1.81502e-05,3.61864e-09,-2.59629e-13,-160343,-46.699], Tmin=(925.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1307.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCsFHH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds)"""),
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
    label = 'CF3(45)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-483.19,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(68.9952,'amu')),
        NonlinearRotor(inertia=([46.5949,46.5949,90.0934],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([504.804,504.824,700.362,1092.66,1279.94,1279.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0058,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1006.05,'J/mol'), sigma=(4.32,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05587,-0.0054938,7.30062e-05,-1.34908e-07,7.857e-11,-58113,8.13986], Tmin=(10,'K'), Tmax=(570.231,'K')), NASAPolynomial(coeffs=[3.53924,0.0118582,-8.75006e-06,2.89362e-09,-3.54041e-13,-58277.3,8.38507], Tmin=(570.231,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=CC(F)[C]=C(F)F(3547)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {6,D} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-460.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.690419,'amu*angstrom^2'), symmetry=1, barrier=(15.8741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690249,'amu*angstrom^2'), symmetry=1, barrier=(15.8702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.584815,0.0832171,-0.000133011,1.18715e-07,-4.22755e-11,-55277,27.2576], Tmin=(100,'K'), Tmax=(793.211,'K')), NASAPolynomial(coeffs=[8.20698,0.0334102,-1.73227e-05,3.41254e-09,-2.39613e-13,-56128.5,-5.49449], Tmin=(793.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-460.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cds_S)"""),
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
    label = 'FC=C(C(F)[C]=C(F)F)C(F)(F)F(12559)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {5,S} {10,D} {15,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {8,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1142.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,219,296,586,564,718,793,1177,1228,350,440,435,1725,194,682,905,1196,1383,3221,562,600,623,1070,1265,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.91477,'amu*angstrom^2'), symmetry=1, barrier=(21.0324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26296,'amu*angstrom^2'), symmetry=1, barrier=(6.04596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262786,'amu*angstrom^2'), symmetry=1, barrier=(6.04197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (207.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36946,0.130986,-0.000223284,1.98926e-07,-6.93143e-11,-137202,37.3961], Tmin=(100,'K'), Tmax=(810.987,'K')), NASAPolynomial(coeffs=[12.8263,0.0416669,-2.23797e-05,4.42636e-09,-3.09866e-13,-138869,-24.203], Tmin=(810.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1142.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(Cds_S)"""),
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
    label = 'F[CH]C(C=C=C(F)F)C(F)(F)F(12560)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
9  C u0 p0 c0 {7,S} {12,D} {14,S}
10 C u1 p0 c0 {4,S} {7,S} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {9,D} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-927.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,94,120,354,641,825,1294,540,610,2055,357.556,357.559,2731.22],'cm^-1')),
        HinderedRotor(inertia=(0.121153,'amu*angstrom^2'), symmetry=1, barrier=(10.9914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121153,'amu*angstrom^2'), symmetry=1, barrier=(10.9915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369095,'amu*angstrom^2'), symmetry=1, barrier=(33.4853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.755503,0.112498,-0.000159729,1.19051e-07,-3.55197e-11,-111359,36.6318], Tmin=(100,'K'), Tmax=(818.822,'K')), NASAPolynomial(coeffs=[15.1004,0.035043,-1.78427e-05,3.53342e-09,-2.50981e-13,-113955,-36.6946], Tmin=(818.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-927.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CsCsF1sH)"""),
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
    label = 'F[CH]C(C(F)=C=C(F)F)C(F)(F)F(12561)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u0 p0 c0 {4,S} {8,S} {13,D}
11 C u1 p0 c0 {5,S} {8,S} {15,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u0 p0 c0 {10,D} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1112.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,145,326,398,834,1303,334,575,1197,1424,3202,94,120,354,641,825,1294,540,610,2055,246.429,248.062,913.928],'cm^-1')),
        HinderedRotor(inertia=(0.0868896,'amu*angstrom^2'), symmetry=1, barrier=(3.73359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00277896,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00279707,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (207.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16741,0.122326,-0.00018211,1.40075e-07,-4.28235e-11,-133570,39.3334], Tmin=(100,'K'), Tmax=(801.702,'K')), NASAPolynomial(coeffs=[16.4159,0.0345934,-1.79539e-05,3.56408e-09,-2.52779e-13,-136390,-41.6089], Tmin=(801.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1112.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCsFHH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CsCsF1sH)"""),
)

species(
    label = 'FC#CC(F)C([CH]F)C(F)(F)F(12562)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
10 C u1 p0 c0 {5,S} {7,S} {15,S}
11 C u0 p0 c0 {8,S} {12,T}
12 C u0 p0 c0 {6,S} {11,T}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-860.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,233,378,609,1068,1270,1314,3037,193,295,551,588,656,1146,1192,1350,334,575,1197,1424,3202,2175,525,239,401,1367,284.094,284.095,2279.87,2279.87],'cm^-1')),
        HinderedRotor(inertia=(0.146155,'amu*angstrom^2'), symmetry=1, barrier=(8.37079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146155,'amu*angstrom^2'), symmetry=1, barrier=(8.37079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682898,'amu*angstrom^2'), symmetry=1, barrier=(39.112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682902,'amu*angstrom^2'), symmetry=1, barrier=(39.112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04501,0.122889,-0.000205429,1.82053e-07,-6.32075e-11,-103300,36.0132], Tmin=(100,'K'), Tmax=(817.123,'K')), NASAPolynomial(coeffs=[11.8405,0.0412554,-2.15103e-05,4.20841e-09,-2.9297e-13,-104787,-19.7591], Tmin=(817.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-860.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    label = 'FC=C([CH][C]=C(F)F)C(F)(F)F(3749)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {8,S} {12,S} {13,S}
10 C u0 p0 c0 {4,S} {8,D} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-838.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,350,440,435,1725,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,562,600,623,1070,1265,1685,370,356.047,356.051,356.051,356.076],'cm^-1')),
        HinderedRotor(inertia=(0.00132977,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132969,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132966,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.241969,0.095507,-0.000113303,6.85108e-08,-1.64715e-11,-100638,34.0453], Tmin=(100,'K'), Tmax=(1012.06,'K')), NASAPolynomial(coeffs=[16.8481,0.0279622,-1.31942e-05,2.56774e-09,-1.82338e-13,-104097,-48.6094], Tmin=(1012.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-838.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)[C]=C(F)F(3068)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {8,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-730.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,234,589,736,816,1240,3237,182,240,577,636,1210,1413,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.500139,'amu*angstrom^2'), symmetry=1, barrier=(11.4992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116281,'amu*angstrom^2'), symmetry=1, barrier=(2.67352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12112,'amu*angstrom^2'), symmetry=1, barrier=(48.7688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3427.25,'J/mol'), sigma=(5.20605,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.33 K, Pc=55.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.866588,0.119074,-0.000202178,1.8159e-07,-6.3977e-11,-87680,37.0298], Tmin=(100,'K'), Tmax=(806.144,'K')), NASAPolynomial(coeffs=[11.336,0.0403044,-2.17058e-05,4.30424e-09,-3.02043e-13,-89055.3,-15.5385], Tmin=(806.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-730.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'FC[C](C(F)[C]=C(F)F)C(F)(F)F(12563)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
10 C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {8,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1063.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20346,0.131057,-0.000231174,2.1892e-07,-8.03645e-11,-127764,41.3728], Tmin=(100,'K'), Tmax=(815.916,'K')), NASAPolynomial(coeffs=[8.29518,0.0526978,-2.86678e-05,5.70026e-09,-4.00116e-13,-128255,3.96476], Tmin=(815.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1063.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Tertalkyl) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH]C([C](F)C=C(F)F)C(F)(F)F(12564)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u1 p0 c0 {4,S} {8,S} {11,S}
11 C u0 p0 c0 {10,S} {13,D} {15,S}
12 C u1 p0 c0 {5,S} {8,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1181.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28772,0.126634,-0.000190014,1.52564e-07,-4.93953e-11,-141904,39.0041], Tmin=(100,'K'), Tmax=(754.137,'K')), NASAPolynomial(coeffs=[14.5929,0.0424127,-2.2516e-05,4.51231e-09,-3.21667e-13,-144300,-33.1316], Tmin=(754.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1181.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCdCsF1s) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'FCC([C](F)[C]=C(F)F)C(F)(F)F(12565)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
10 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
11 C u1 p0 c0 {5,S} {8,S} {13,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {11,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
"""),
    E0 = (-1111.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.79055,0.115605,-0.000155753,1.13121e-07,-3.36323e-11,-133481,38.9746], Tmin=(100,'K'), Tmax=(813.697,'K')), NASAPolynomial(coeffs=[13.7999,0.0438788,-2.35264e-05,4.78352e-09,-3.45733e-13,-135856,-28.4075], Tmin=(813.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1111.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCdCsF1s) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47624,0.139156,-0.000253687,2.40575e-07,-8.71736e-11,-136196,40.6148], Tmin=(100,'K'), Tmax=(835.539,'K')), NASAPolynomial(coeffs=[8.73412,0.0518702,-2.8041e-05,5.52236e-09,-3.8396e-13,-136562,1.21255], Tmin=(835.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1133.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[C](F)[C]=CC(C(F)F)C(F)(F)F(7283)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u0 p0 c0 {8,S} {13,D} {16,S}
12 C u1 p0 c0 {6,S} {7,S} {13,S}
13 C u1 p0 c0 {11,D} {12,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1127.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884241,0.114225,-0.000143266,9.22613e-08,-2.38398e-11,-135421,39.6903], Tmin=(100,'K'), Tmax=(938.904,'K')), NASAPolynomial(coeffs=[17.4323,0.0361911,-1.85986e-05,3.74182e-09,-2.69899e-13,-138861,-47.5216], Tmin=(938.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1127.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsH))"""),
)

species(
    label = 'F[CH]C([CH]C(F)=C(F)F)C(F)(F)F(12567)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u1 p0 c0 {8,S} {11,S} {15,S}
11 C u0 p0 c0 {4,S} {10,S} {13,D}
12 C u1 p0 c0 {5,S} {8,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1147.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30035,0.126496,-0.000191408,1.54987e-07,-5.04756e-11,-137793,37.7921], Tmin=(100,'K'), Tmax=(750.847,'K')), NASAPolynomial(coeffs=[14.6415,0.0415636,-2.17244e-05,4.318e-09,-3.06068e-13,-140187,-34.5491], Tmin=(750.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1147.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(CsCsFHH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Allyl_S) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[C](F)C(C(F)F)C(F)[C]=C(F)F(6744)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {13,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {16,S}
11 C u1 p0 c0 {4,S} {5,S} {8,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1024.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,562,600,623,1070,1265,1685,370,198.829,198.908,198.993,1741.43],'cm^-1')),
        HinderedRotor(inertia=(0.174321,'amu*angstrom^2'), symmetry=1, barrier=(4.89806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174391,'amu*angstrom^2'), symmetry=1, barrier=(4.89899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1745,'amu*angstrom^2'), symmetry=1, barrier=(4.8984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16607,'amu*angstrom^2'), symmetry=1, barrier=(32.7679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3374.93,'J/mol'), sigma=(5.29414,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.16 K, Pc=51.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30365,0.128321,-0.00020668,1.81642e-07,-6.41986e-11,-123062,42.2379], Tmin=(100,'K'), Tmax=(757.094,'K')), NASAPolynomial(coeffs=[12.4038,0.046301,-2.51612e-05,5.05854e-09,-3.59765e-13,-124863,-18.2616], Tmin=(757.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1024.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH]C([C](F)F)C(F)C(F)=C(F)F(7286)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {11,S} {12,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
10 C u0 p0 c0 {2,S} {9,S} {13,D}
11 C u1 p0 c0 {3,S} {8,S} {16,S}
12 C u1 p0 c0 {4,S} {5,S} {8,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1032.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,323,467,575,827,1418,334,575,1197,1424,3202,190,488,555,1236,1407,182,240,577,636,1210,1413,183.668,184.881,191.168,2141.17],'cm^-1')),
        HinderedRotor(inertia=(0.214089,'amu*angstrom^2'), symmetry=1, barrier=(5.08865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220148,'amu*angstrom^2'), symmetry=1, barrier=(5.23267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215416,'amu*angstrom^2'), symmetry=1, barrier=(5.1205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71718,'amu*angstrom^2'), symmetry=1, barrier=(40.8865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80046,0.142647,-0.000248344,2.23902e-07,-7.82901e-11,-123946,41.5549], Tmin=(100,'K'), Tmax=(824.211,'K')), NASAPolynomial(coeffs=[12.9591,0.0452758,-2.42892e-05,4.78298e-09,-3.33232e-13,-125505,-21.494], Tmin=(824.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1032.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFHH) + group(CsCsFFH) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]=[C]C(F)C(C(F)F)C(F)(F)F(7287)',
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
10 C u0 p0 c0 {2,S} {3,S} {8,S} {16,S}
11 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
12 C u1 p0 c0 {9,S} {13,D}
13 C u1 p0 c0 {7,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1016.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,235,523,627,1123,1142,1372,1406,3097,193,295,551,588,656,1146,1192,1350,1685,370,167,640,1190,281.441,281.634,2215.71],'cm^-1')),
        HinderedRotor(inertia=(0.143366,'amu*angstrom^2'), symmetry=1, barrier=(8.05834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143216,'amu*angstrom^2'), symmetry=1, barrier=(8.05662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143345,'amu*angstrom^2'), symmetry=1, barrier=(8.05813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.596047,'amu*angstrom^2'), symmetry=1, barrier=(33.5304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01004,0.124375,-0.000175772,1.14466e-07,-1.98418e-11,-122140,38.7855], Tmin=(100,'K'), Tmax=(586.235,'K')), NASAPolynomial(coeffs=[14.2217,0.042495,-2.26845e-05,4.53309e-09,-3.21677e-13,-124305,-29.7962], Tmin=(586.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1016.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCsFFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]=C(F)C(F)C([CH]F)C(F)(F)F(12568)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
11 C u0 p0 c0 {5,S} {9,S} {13,D}
12 C u1 p0 c0 {6,S} {8,S} {16,S}
13 C u1 p0 c0 {7,S} {11,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1030.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,193,295,551,588,656,1146,1192,1350,246,474,533,1155,334,575,1197,1424,3202,167,640,1190,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69562,0.138483,-0.000233996,2.07009e-07,-7.17503e-11,-123798,40.4679], Tmin=(100,'K'), Tmax=(809.963,'K')), NASAPolynomial(coeffs=[13.6298,0.043765,-2.33355e-05,4.60517e-09,-3.22117e-13,-125656,-26.3842], Tmin=(809.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1030.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFF) + group(CsCsFHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Cdj(Cd-CsF1s)(F1s))"""),
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
    E0 = (-380.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (157.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-210.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-205.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-286.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (57.7147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (80.7126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-372.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-291.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-302.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-305.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-249.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-254.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-120.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-260.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-224.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-71.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-84.2179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-147.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-157.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-18.4574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-212.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-197.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-224.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-238.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-187.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-240.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-136.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-126.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-162.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-186.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['FC=C=C(F)F(1325)', 'FC=CC(F)(F)F(1027)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'F[CH]C(F)C(F)[C]=C(F)F(3030)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33582e-06,'m^3/(mol*s)'), n=3.3552, Ea=(235.61,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R',), comment="""Estimated from node CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC(F)=[C]C(F)C(F)[CH]C(F)(F)F(12466)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['FC(F)=[C]C(F)[CH]C(F)C(F)(F)F(12555)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['F[CH]C(=C(F)F)C([CH]F)C(F)(F)F(7322)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(804)', 'FC(F)=[C]C(F)[CH]C(F)(F)F(3367)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=C(F)F(1205)', 'F[CH]C([CH]F)C(F)(F)F(12266)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['FC(F)=C1C(F)C(C(F)(F)F)C1F(12556)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_1H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['FC=C(C(F)C=C(F)F)C(F)(F)F(12557)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['FCC(C(F)=C=C(F)F)C(F)(F)F(12558)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH][C]=C(F)F(2842)', 'FC=CC(F)(F)F(1027)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(51.3428,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CF3(45)', 'FC=CC(F)[C]=C(F)F(3547)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(18.1875,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'FC=C(C(F)[C]=C(F)F)C(F)(F)F(12559)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'F[CH]C(C=C=C(F)F)C(F)(F)F(12560)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(58.1322,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FC=C=C(F)F(1325)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(9.20217,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'F[CH]C(C(F)=C=C(F)F)C(F)(F)F(12561)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'FC#CC(F)C([CH]F)C(F)(F)F(12562)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(21.8547,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'FC=C([CH][C]=C(F)F)C(F)(F)F(3749)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(295.637,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'F[CH]C(=C(F)F)C(F)[C]=C(F)F(3068)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(26.4943,'m^3/(mol*s)'), n=1.22463, Ea=(177.77,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CHF(40)', 'FC(F)=[C]C(F)[CH]C(F)(F)F(3367)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['FC[C](C(F)[C]=C(F)F)C(F)(F)F(12563)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.06147e+10,'s^-1'), n=0.76, Ea=(168.197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['F[CH]C([C](F)C=C(F)F)C(F)(F)F(12564)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['FCC([C](F)[C]=C(F)F)C(F)(F)F(12565)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(91.367,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['F[CH][C](C(F)C=C(F)F)C(F)(F)F(12566)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['F[C](F)[C]=CC(C(F)F)C(F)(F)F(7283)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(193.415,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    products = ['F[CH]C([CH]C(F)=C(F)F)C(F)(F)F(12567)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(140.299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[C](F)C(C(F)F)C(F)[C]=C(F)F(6744)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(212.004,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[CH]C([C](F)F)C(F)C(F)=C(F)F(7286)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(229.457,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[C]=[C]C(F)C(C(F)F)C(F)(F)F(7287)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(178.726,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]=C(F)C(F)C([CH]F)C(F)(F)F(12568)'],
    products = ['F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3689',
    isomers = [
        'F[CH]C(C(F)[C]=C(F)F)C(F)(F)F(7285)',
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
    label = 'PDepNetwork #3689',
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

