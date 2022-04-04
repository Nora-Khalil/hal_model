species(
    label = 'F[CH]C(C[C](F)F)C(F)(F)F(6814)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
10 C u1 p0 c0 {4,S} {7,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1036.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,193,295,551,588,656,1146,1192,1350,334,575,1197,1424,3202,190,488,555,1236,1407,281.01,281.052,1818.06,1818.09],'cm^-1')),
        HinderedRotor(inertia=(0.168025,'amu*angstrom^2'), symmetry=1, barrier=(9.41948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168029,'amu*angstrom^2'), symmetry=1, barrier=(9.4197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.641767,'amu*angstrom^2'), symmetry=1, barrier=(35.9691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16802,'amu*angstrom^2'), symmetry=1, barrier=(9.41818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.703907,0.113643,-0.000182112,1.59065e-07,-5.53581e-11,-124540,35.5354], Tmin=(100,'K'), Tmax=(790.204,'K')), NASAPolynomial(coeffs=[11.3073,0.0408724,-2.12553e-05,4.18659e-09,-2.93824e-13,-126064,-17.2187], Tmin=(790.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1036.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFFF) + group(CsCsFHH) + group(CsCsFFH) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'CH2CF2(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.616,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[CH]C(F)C[C](F)F(802)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {5,S} {12,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-553.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.217836,'amu*angstrom^2'), symmetry=1, barrier=(5.00848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217911,'amu*angstrom^2'), symmetry=1, barrier=(5.01019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547638,'amu*angstrom^2'), symmetry=1, barrier=(12.5913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2936.11,'J/mol'), sigma=(5.04912,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.61 K, Pc=51.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288231,0.0947533,-0.000173412,1.65289e-07,-5.99531e-11,-66460.3,27.1491], Tmin=(100,'K'), Tmax=(839.561,'K')), NASAPolynomial(coeffs=[6.92912,0.0356971,-1.89165e-05,3.71359e-09,-2.57782e-13,-66609.1,2.02637], Tmin=(839.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-553.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)CC(F)[CH]C(F)(F)F(12102)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1045.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.581296,0.111401,-0.000179218,1.59673e-07,-5.69877e-11,-125537,34.65], Tmin=(100,'K'), Tmax=(777.116,'K')), NASAPolynomial(coeffs=[10.2117,0.0431422,-2.29421e-05,4.57004e-09,-3.2312e-13,-126830,-12.2302], Tmin=(777.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1045.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFFF) + group(CsCsFFH) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)C[CH]C(F)C(F)(F)F(12132)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {10,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
10 C u1 p0 c0 {7,S} {8,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1025.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,193,295,551,588,656,1146,1192,1350,3025,407.5,1350,352.5,190,488,555,1236,1407,242.704,242.764,1238.02,1238.05],'cm^-1')),
        HinderedRotor(inertia=(0.15456,'amu*angstrom^2'), symmetry=1, barrier=(6.46578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154729,'amu*angstrom^2'), symmetry=1, barrier=(6.46546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154659,'amu*angstrom^2'), symmetry=1, barrier=(6.46576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154576,'amu*angstrom^2'), symmetry=1, barrier=(6.46577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.665723,0.113579,-0.000185609,1.6508e-07,-5.81735e-11,-123133,36.4237], Tmin=(100,'K'), Tmax=(798.482,'K')), NASAPolynomial(coeffs=[10.6039,0.0417765,-2.18917e-05,4.31796e-09,-3.0292e-13,-124444,-12.3456], Tmin=(798.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1025.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCsFFH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(F)C([CH]F)C(F)(F)F(6839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
10 C u1 p0 c0 {6,S} {7,S} {13,S}
11 C u1 p0 c0 {8,S} {14,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1050.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.850315,0.115373,-0.000171876,1.3517e-07,-4.25688e-11,-126159,34.2438], Tmin=(100,'K'), Tmax=(776.827,'K')), NASAPolynomial(coeffs=[14.5573,0.036039,-1.8691e-05,3.71078e-09,-2.63149e-13,-128553,-36.1983], Tmin=(776.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1050.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCsFHH) + group(Cs-CsHHH) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = 'F[C](F)C[CH]C(F)(F)F(1048)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-827.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,1704.75],'cm^-1')),
        HinderedRotor(inertia=(0.162916,'amu*angstrom^2'), symmetry=1, barrier=(3.74575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163006,'amu*angstrom^2'), symmetry=1, barrier=(3.74783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1631,'amu*angstrom^2'), symmetry=1, barrier=(3.74999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615058,0.0827085,-0.000133047,1.19867e-07,-4.30679e-11,-99445.5,30.5336], Tmin=(100,'K'), Tmax=(787.654,'K')), NASAPolynomial(coeffs=[8.05947,0.0334273,-1.73425e-05,3.43745e-09,-2.4245e-13,-100262,-1.34474], Tmin=(787.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-827.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = '[CH2]C([CH]F)C(F)(F)F(12133)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
7  C u1 p0 c0 {4,S} {5,S} {10,S}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-602.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,180,1247.35],'cm^-1')),
        HinderedRotor(inertia=(0.101819,'amu*angstrom^2'), symmetry=1, barrier=(2.34102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101309,'amu*angstrom^2'), symmetry=1, barrier=(2.32929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102395,'amu*angstrom^2'), symmetry=1, barrier=(2.35427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.947415,0.0715373,-9.63716e-05,7.27073e-08,-2.22302e-11,-72366.3,25.871], Tmin=(100,'K'), Tmax=(797.486,'K')), NASAPolynomial(coeffs=[9.87259,0.0267687,-1.21623e-05,2.30852e-09,-1.60234e-13,-73789.8,-15.1677], Tmin=(797.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-602.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Isobutyl)"""),
)

species(
    label = 'FC1C(C(F)(F)F)CC1(F)F(12134)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {7,S} {10,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
11 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-1304.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.027077,0.0795657,-6.79941e-05,2.84197e-08,-4.70443e-12,-156802,29.5728], Tmin=(100,'K'), Tmax=(1445.06,'K')), NASAPolynomial(coeffs=[19.8629,0.0245091,-1.08442e-05,2.05397e-09,-1.43078e-13,-162551,-73.7074], Tmin=(1445.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1304.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsFFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FC=C(CC(F)F)C(F)(F)F(12135)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {14,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1300.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.155005,0.0978194,-0.000121363,8.0041e-08,-2.14325e-11,-156253,30.7188], Tmin=(100,'K'), Tmax=(903.269,'K')), NASAPolynomial(coeffs=[13.9456,0.0353775,-1.76715e-05,3.51102e-09,-2.51347e-13,-158800,-35.8742], Tmin=(903.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1300.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(CdCFH)"""),
)

species(
    label = 'FCC(C=C(F)F)C(F)(F)F(12136)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
9  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
10 C u0 p0 c0 {7,S} {11,D} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1287.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0281126,0.0914621,-0.000102232,5.96475e-08,-1.41154e-11,-154678,31.4851], Tmin=(100,'K'), Tmax=(1016.14,'K')), NASAPolynomial(coeffs=[14.8326,0.0331854,-1.62068e-05,3.2087e-09,-2.29937e-13,-157687,-40.1753], Tmin=(1016.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1287.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF)"""),
)

species(
    label = '[CH2][C](F)F(141)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25818,0.0186259,-1.56413e-05,8.22305e-09,-2.04579e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.66,'K')), NASAPolynomial(coeffs=[4.47009,0.0131647,-6.4128e-06,1.29206e-09,-9.37476e-14,-13305.9,7.85289], Tmin=(887.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
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
    label = 'FC=CC[C](F)F(809)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {7,D} {10,S}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u0 p0 c0 {3,S} {5,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-426.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,194,682,905,1196,1383,3221,180,1013.33],'cm^-1')),
        HinderedRotor(inertia=(0.23822,'amu*angstrom^2'), symmetry=1, barrier=(5.47714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236881,'amu*angstrom^2'), symmetry=1, barrier=(5.44637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83808,0.0381015,-2.41359e-05,6.97405e-09,-7.1109e-13,-51325.6,13.629], Tmin=(10,'K'), Tmax=(1627.49,'K')), NASAPolynomial(coeffs=[18.1768,0.00841183,-1.88875e-06,-4.30274e-11,4.48415e-14,-56728.1,-64.7895], Tmin=(1627.49,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""FCDCC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC=C(C[C](F)F)C(F)(F)F(12137)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {11,D}
10 C u1 p0 c0 {4,S} {5,S} {7,S}
11 C u0 p0 c0 {6,S} {9,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-1099.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,219,296,586,564,718,793,1177,1228,350,440,435,1725,190,488,555,1236,1407,194,682,905,1196,1383,3221,395.391,396.664,2652.48],'cm^-1')),
        HinderedRotor(inertia=(0.0805752,'amu*angstrom^2'), symmetry=1, barrier=(9.14137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291487,'amu*angstrom^2'), symmetry=1, barrier=(32.1979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0824025,'amu*angstrom^2'), symmetry=1, barrier=(9.14483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00414646,0.0940735,-0.000117535,7.67088e-08,-2.02104e-11,-132051,32.1099], Tmin=(100,'K'), Tmax=(919.319,'K')), NASAPolynomial(coeffs=[14.265,0.0319879,-1.62341e-05,3.24803e-09,-2.33482e-13,-134674,-35.5303], Tmin=(919.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1099.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
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
    label = 'F[CH]C(C=C(F)F)C(F)(F)F(12138)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {4,S} {7,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {9,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-1091.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,193,295,551,588,656,1146,1192,1350,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,182,240,577,636,1210,1413,353.595,353.727,2864.68],'cm^-1')),
        HinderedRotor(inertia=(0.125321,'amu*angstrom^2'), symmetry=1, barrier=(11.1328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125394,'amu*angstrom^2'), symmetry=1, barrier=(11.1373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.4192,'amu*angstrom^2'), symmetry=1, barrier=(37.2463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.226492,0.0999259,-0.000137051,9.9126e-08,-2.88132e-11,-131189,33.0817], Tmin=(100,'K'), Tmax=(838.477,'K')), NASAPolynomial(coeffs=[13.8993,0.0325374,-1.64949e-05,3.27179e-09,-2.33053e-13,-133558,-32.5785], Tmin=(838.477,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1091.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFF) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCsF1sH)"""),
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
    label = 'F[CH]C(C[C](F)F)=C(F)F(3001)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {9,S} {10,D}
8  C u1 p0 c0 {1,S} {2,S} {6,S}
9  C u1 p0 c0 {3,S} {7,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {7,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-687.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,190,488,555,1236,1407,234,589,736,816,1240,3237,182,240,577,636,1210,1413,476.789,476.933],'cm^-1')),
        HinderedRotor(inertia=(0.226691,'amu*angstrom^2'), symmetry=1, barrier=(36.5718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00113527,'amu*angstrom^2'), symmetry=1, barrier=(7.72525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226546,'amu*angstrom^2'), symmetry=1, barrier=(36.5726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3184,'J/mol'), sigma=(5.13001,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=497.33 K, Pc=53.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5064,0.0821003,-9.63392e-05,5.9386e-08,-1.49102e-11,-82529.6,31.714], Tmin=(100,'K'), Tmax=(958.189,'K')), NASAPolynomial(coeffs=[12.8858,0.0304206,-1.54352e-05,3.09515e-09,-2.23028e-13,-84901.9,-27.4806], Tmin=(958.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cs-CdHH)(F1s)(F1s)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'FC[C](C[C](F)F)C(F)(F)F(12139)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {10,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
9  C u0 p0 c0 {4,S} {10,S} {14,S} {15,S}
10 C u1 p0 c0 {7,S} {8,S} {9,S}
11 C u1 p0 c0 {5,S} {6,S} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-1044.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.689432,0.118472,-0.000209777,1.97102e-07,-7.10867e-11,-125445,35.8927], Tmin=(100,'K'), Tmax=(837.552,'K')), NASAPolynomial(coeffs=[7.8783,0.0462093,-2.42247e-05,4.72557e-09,-3.27227e-13,-125781,2.64032], Tmin=(837.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1044.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFFF) + group(CsCsFHH) + group(CsCsFFH) + radical(Tertalkyl) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C([CH]C(F)F)C(F)(F)F(12140)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u1 p0 c0 {6,S} {7,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1042.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.541802,0.10938,-0.000170179,1.47124e-07,-5.14328e-11,-125267,36.1495], Tmin=(100,'K'), Tmax=(759.071,'K')), NASAPolynomial(coeffs=[10.9937,0.0414284,-2.17412e-05,4.32174e-09,-3.05821e-13,-126811,-14.9634], Tmin=(759.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1042.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFFF) + group(CsCsFHH) + group(CsCsFFH) + radical(Cs_S) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'FCC([CH][C](F)F)C(F)(F)F(12141)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
9  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
10 C u1 p0 c0 {7,S} {11,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1035.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.512706,0.108953,-0.000171148,1.4933e-07,-5.24571e-11,-124349,36.6826], Tmin=(100,'K'), Tmax=(771.388,'K')), NASAPolynomial(coeffs=[10.6545,0.0415838,-2.17472e-05,4.31062e-09,-3.04263e-13,-125790,-12.4695], Tmin=(771.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1035.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFFF) + group(CsCsFHH) + group(CsCsFFH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][C](CC(F)F)C(F)(F)F(12142)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {14,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {11,S}
11 C u1 p0 c0 {6,S} {10,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1051.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.733318,0.11909,-0.00020955,1.95956e-07,-7.055e-11,-126362,35.4116], Tmin=(100,'K'), Tmax=(834.133,'K')), NASAPolynomial(coeffs=[8.26744,0.0459634,-2.41641e-05,4.72337e-09,-3.27652e-13,-126821,-0.130886], Tmin=(834.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1051.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFFF) + group(CsCsFHH) + group(CsCsFFH) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[C](F)CC([C](F)F)C(F)F(6677)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {15,S}
10 C u1 p0 c0 {3,S} {4,S} {7,S}
11 C u1 p0 c0 {5,S} {6,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-1005.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,146,234,414,562,504,606,1176,1296,1354,1460,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.100989,'amu*angstrom^2'), symmetry=1, barrier=(2.32193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103623,'amu*angstrom^2'), symmetry=1, barrier=(2.38249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10233,'amu*angstrom^2'), symmetry=1, barrier=(2.35277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926774,'amu*angstrom^2'), symmetry=1, barrier=(21.3084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3132.16,'J/mol'), sigma=(5.21338,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.24 K, Pc=50.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.822414,0.116146,-0.000186804,1.61828e-07,-5.57295e-11,-120742,36.8737], Tmin=(100,'K'), Tmax=(791.505,'K')), NASAPolynomial(coeffs=[12.1324,0.0395501,-2.05604e-05,4.04544e-09,-2.83618e-13,-122444,-20.3951], Tmin=(791.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1005.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFFH) + group(CsCsFFH) + group(CsCsFFH) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(CC(F)(F)F)[C](F)F(6815)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u1 p0 c0 {4,S} {7,S} {15,S}
11 C u1 p0 c0 {5,S} {6,S} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-1038.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (178.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22927,0.104887,-0.000140804,8.68073e-08,-1.22849e-11,-124745,34.0083], Tmin=(100,'K'), Tmax=(584.234,'K')), NASAPolynomial(coeffs=[12.0875,0.0392481,-2.026e-05,4.00779e-09,-2.83302e-13,-126503,-21.523], Tmin=(584.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1038.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFFF) + radical(Csj(Cs-CsCsH)(F1s)(H)) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
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
    E0 = (-402.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (121.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-242.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-230.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-156.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (21.3978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (65.4137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-394.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-339.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-339.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-285.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-260.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-253.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-295.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-219.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-55.9212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-166.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-54.7743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-169.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-234.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-238.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-246.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-280.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-159.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-161.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['CH2CF2(57)', 'FC=CC(F)(F)F(1027)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'F[CH]C(F)C[C](F)F(802)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33582e-06,'m^3/(mol*s)'), n=3.3552, Ea=(244.343,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R',), comment="""Estimated from node CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['F[C](F)CC(F)[CH]C(F)(F)F(12102)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)C[CH]C(F)C(F)(F)F(12132)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['[CH2]C(F)(F)C([CH]F)C(F)(F)F(6839)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(804)', 'F[C](F)C[CH]C(F)(F)F(1048)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]F(138)', '[CH2]C([CH]F)C(F)(F)F(12133)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['FC1C(C(F)(F)F)CC1(F)F(12134)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_1H;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['FC=C(CC(F)F)C(F)(F)F(12135)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['FCC(C=C(F)F)C(F)(F)F(12136)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C](F)F(141)', 'FC=CC(F)(F)F(1027)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(21.1625,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CF3(45)', 'FC=CC[C](F)F(809)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(15.1465,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'FC=C(C[C](F)F)C(F)(F)F(12137)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2CF2(57)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(13.1888,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'F[CH]C(C=C(F)F)C(F)(F)F(12138)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(306600,'m^3/(mol*s)'), n=0.481, Ea=(26.1491,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C](F)F(141)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -1.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'F[CH]C(C[C](F)F)=C(F)F(3001)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(26.4943,'m^3/(mol*s)'), n=1.22463, Ea=(168.042,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHF(40)', 'F[C](F)C[CH]C(F)(F)F(1048)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CF2(43)', '[CH2]C([CH]F)C(F)(F)F(12133)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(2.25929,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['FC[C](C[C](F)F)C(F)(F)F(12139)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.06147e+10,'s^-1'), n=0.76, Ea=(168.197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['F[CH]C([CH]C(F)F)C(F)(F)F(12140)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.23544e+08,'s^-1'), n=1.37575, Ea=(164.097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] + [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['FCC([CH][C](F)F)C(F)(F)F(12141)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(182.734,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['F[CH][C](CC(F)F)C(F)(F)F(12142)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.50974e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C](F)CC([C](F)F)C(F)F(6677)'],
    products = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(212.004,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[CH]C(C[C](F)F)C(F)(F)F(6814)'],
    products = ['F[CH]C(CC(F)(F)F)[C](F)F(6815)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0108995,'s^-1'), n=4.43046, Ea=(241.482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #3633',
    isomers = [
        'F[CH]C(C[C](F)F)C(F)(F)F(6814)',
    ],
    reactants = [
        ('CH2CF2(57)', 'FC=CC(F)(F)F(1027)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3633',
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

