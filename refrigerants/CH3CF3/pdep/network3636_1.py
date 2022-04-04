species(
    label = '[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1151.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,351,323,533,609,664,892,1120,1201,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,180,180,1187.18],'cm^-1')),
        HinderedRotor(inertia=(0.262837,'amu*angstrom^2'), symmetry=1, barrier=(6.04314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265249,'amu*angstrom^2'), symmetry=1, barrier=(6.09861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264873,'amu*angstrom^2'), symmetry=1, barrier=(6.08996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.192276,0.100675,-0.000150699,1.19286e-07,-3.80135e-11,-138394,32.0787], Tmin=(100,'K'), Tmax=(766.062,'K')), NASAPolynomial(coeffs=[12.8104,0.0327806,-1.77577e-05,3.59345e-09,-2.57764e-13,-140386,-27.1867], Tmin=(766.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1151.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H))"""),
)

species(
    label = 'CF2O(49)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-618.61,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C([CH]F)C(F)(F)F(6991)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {8,S} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1154.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,193,295,551,588,656,1146,1192,1350,334,575,1197,1424,3202,264.244,264.783,2653.69],'cm^-1')),
        HinderedRotor(inertia=(0.147282,'amu*angstrom^2'), symmetry=1, barrier=(7.38007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150624,'amu*angstrom^2'), symmetry=1, barrier=(7.38536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.670751,'amu*angstrom^2'), symmetry=1, barrier=(33.5695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3610.62,'J/mol'), sigma=(6.35752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.97 K, Pc=31.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.382049,0.104333,-0.000158315,1.24894e-07,-3.92798e-11,-138721,32.6969], Tmin=(100,'K'), Tmax=(778.899,'K')), NASAPolynomial(coeffs=[13.9897,0.0305255,-1.6175e-05,3.23338e-09,-2.29991e-13,-140960,-33.0473], Tmin=(778.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1154.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCsFHH) + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
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
    label = '[O]C(F)(F)[CH]F(841)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u1 p0 c0 {3,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-445.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0874832,'amu*angstrom^2'), symmetry=1, barrier=(9.70476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3211.42,'J/mol'), sigma=(5.51205,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.62 K, Pc=43.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99253,0.0473334,-6.85919e-05,5.03377e-08,-1.45679e-11,-53519.7,18.0995], Tmin=(100,'K'), Tmax=(847.739,'K')), NASAPolynomial(coeffs=[9.47107,0.0120471,-6.15723e-06,1.23976e-09,-8.91187e-14,-54787.7,-16.745], Tmin=(847.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
)

species(
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,29230.2,5.12617], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C](F)C(F)[CH]C(F)(F)F(1140)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1004.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,738.616],'cm^-1')),
        HinderedRotor(inertia=(0.203364,'amu*angstrom^2'), symmetry=1, barrier=(4.67574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203368,'amu*angstrom^2'), symmetry=1, barrier=(4.67583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203349,'amu*angstrom^2'), symmetry=1, barrier=(4.6754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0384801,0.0960041,-0.000159619,1.40805e-07,-4.90751e-11,-120719,30.6515], Tmin=(100,'K'), Tmax=(788.794,'K')), NASAPolynomial(coeffs=[10.7976,0.0311133,-1.65743e-05,3.30396e-09,-2.33216e-13,-122094,-16.6653], Tmin=(788.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1004.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'FC1C(C(F)(F)F)OC1(F)F(12107)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
11 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1461.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433412,0.0775935,-7.79553e-05,3.89673e-08,-7.7818e-12,-175605,27.4385], Tmin=(100,'K'), Tmax=(1201.91,'K')), NASAPolynomial(coeffs=[16.3163,0.024734,-1.19852e-05,2.37501e-09,-1.7044e-13,-179423,-52.1083], Tmin=(1201.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1461.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + ring(Cs-Cs(F)(F)-O2s-Cs)"""),
)

species(
    label = 'OC(F)(F)C(F)=CC(F)(F)F(7536)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {8,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
10 C u0 p0 c0 {6,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {12,S}
12 H u0 p0 c0 {11,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1478.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.170436,0.0980647,-0.00013768,9.90768e-08,-2.83344e-11,-177677,30.4079], Tmin=(100,'K'), Tmax=(855.317,'K')), NASAPolynomial(coeffs=[14.8869,0.0276484,-1.41908e-05,2.82571e-09,-2.01629e-13,-180253,-39.8819], Tmin=(855.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1478.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH)"""),
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
    label = '[O]C(F)(F)C=CC(F)(F)F(12255)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
8  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {11,S}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1044.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,526,555,698,907,1200,1145,1227,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.376373,'amu*angstrom^2'), symmetry=1, barrier=(8.65357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809698,'amu*angstrom^2'), symmetry=1, barrier=(18.6166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.669218,0.0750653,-8.24579e-05,4.54629e-08,-1.00303e-11,-125540,27.8938], Tmin=(100,'K'), Tmax=(1092.6,'K')), NASAPolynomial(coeffs=[14.7188,0.0236299,-1.1844e-05,2.37697e-09,-1.71766e-13,-128610,-41.1318], Tmin=(1092.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1044.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CsCdFFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(O2sj(Cs-F1sF1sCd))"""),
)

species(
    label = '[O][C](F)F(263)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-255.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53788e-05,4.83907e-09,3.86802e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58999,0.00108969,-4.53602e-07,1.24872e-10,-1.13713e-14,-32130.4,-16.3888], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
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
    label = '[O]C(F)(F)C(F)=CC(F)(F)F(12256)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
10 C u0 p0 c0 {8,S} {11,D} {12,S}
11 C u0 p0 c0 {6,S} {9,S} {10,D}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1230.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,526,555,698,907,1200,1145,1227,3010,987.5,1337.5,450,1655,323,467,575,827,1418,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.176377,'amu*angstrom^2'), symmetry=1, barrier=(4.05526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833195,'amu*angstrom^2'), symmetry=1, barrier=(19.1568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.20834,0.0848099,-0.000101029,5.94185e-08,-1.37786e-11,-147857,30.5493], Tmin=(100,'K'), Tmax=(1050.88,'K')), NASAPolynomial(coeffs=[16.8227,0.02157,-1.07618e-05,2.154e-09,-1.55542e-13,-151349,-50.4298], Tmin=(1050.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1230.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(O2sj(Cs-F1sF1sCd))"""),
)

species(
    label = 'O=C(F)C(F)[CH]C(F)(F)F(12257)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u0 p0 c0 {5,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1124.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,233,1109,1254,1325,1339,3290,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,486,617,768,1157,1926,321.977,321.985,2325.24],'cm^-1')),
        HinderedRotor(inertia=(0.128407,'amu*angstrom^2'), symmetry=1, barrier=(9.44757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128336,'amu*angstrom^2'), symmetry=1, barrier=(9.4473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.392963,'amu*angstrom^2'), symmetry=1, barrier=(28.8901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661762,0.0800049,-0.000105382,7.40165e-08,-2.11788e-11,-135098,29.8789], Tmin=(100,'K'), Tmax=(845.8,'K')), NASAPolynomial(coeffs=[11.3871,0.0292847,-1.54354e-05,3.12359e-09,-2.25407e-13,-136912,-20.0689], Tmin=(845.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1124.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFF) + group(COCsFO) + radical(CCJCC=O)"""),
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
    label = '[O]C(F)(F)C(F)C=C(F)F(3661)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {7,S} {10,D} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-949.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,351,323,533,609,664,892,1120,1201,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.235448,'amu*angstrom^2'), symmetry=1, barrier=(5.41341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236004,'amu*angstrom^2'), symmetry=1, barrier=(5.42619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.348449,0.0874545,-0.00012609,9.67919e-08,-3.00033e-11,-114117,28.883], Tmin=(100,'K'), Tmax=(786.186,'K')), NASAPolynomial(coeffs=[11.7228,0.0295813,-1.5667e-05,3.15259e-09,-2.25845e-13,-115905,-23.2552], Tmin=(786.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-949.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = '[O]C(F)=C[CH]C(F)(F)F(12258)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {10,S}
8  C u0 p0 c0 {7,S} {9,D} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-785.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180,841.827,1277.95],'cm^-1')),
        HinderedRotor(inertia=(0.295256,'amu*angstrom^2'), symmetry=1, barrier=(6.78853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84154,'amu*angstrom^2'), symmetry=1, barrier=(42.3407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21051,0.0638685,-6.64162e-05,3.52287e-08,-7.59499e-12,-94344.4,24.0073], Tmin=(100,'K'), Tmax=(1105.75,'K')), NASAPolynomial(coeffs=[12.3053,0.0237336,-1.19714e-05,2.40351e-09,-1.73508e-13,-96798,-30.634], Tmin=(1105.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-785.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + group(Cds-CdsCsH) + group(CdCFO) + radical(C=COJ) + radical(Allyl_S)"""),
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
    label = 'O=C(F)[C](F)[CH]C(F)(F)F(12259)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {7,S} {9,S} {11,S}
9  C u1 p0 c0 {4,S} {8,S} {10,S}
10 C u0 p0 c0 {5,S} {6,D} {9,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-967.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,234,347,1316,1464,611,648,830,1210,1753,332.49,333.056,2352.18],'cm^-1')),
        HinderedRotor(inertia=(0.134126,'amu*angstrom^2'), symmetry=1, barrier=(10.5661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133882,'amu*angstrom^2'), symmetry=1, barrier=(10.5613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.510528,'amu*angstrom^2'), symmetry=1, barrier=(40.2404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.760417,0.0768558,-0.000102361,7.08349e-08,-1.97246e-11,-116293,30.021], Tmin=(100,'K'), Tmax=(872.727,'K')), NASAPolynomial(coeffs=[12.0922,0.0249182,-1.30926e-05,2.64354e-09,-1.90532e-13,-118271,-23.1056], Tmin=(872.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-967.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFF) + group(COCsFO) + radical(CCJCC=O) + radical(CsCOCsF1s)"""),
)

species(
    label = '[O]C(F)(F)[C](F)CC(F)(F)F(12260)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {9,S} {11,S} {12,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
11 C u1 p0 c0 {6,S} {8,S} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-1163.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.584139,0.110984,-0.000187655,1.64451e-07,-5.64124e-11,-139817,34.4381], Tmin=(100,'K'), Tmax=(809.899,'K')), NASAPolynomial(coeffs=[12.4347,0.03292,-1.75805e-05,3.4689e-09,-2.42517e-13,-141474,-22.838], Tmin=(809.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1163.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + radical(O2sj(Cs-CsF1sF1s)) + radical(CsCsCsF1s)"""),
)

species(
    label = 'OC(F)(F)[C](F)[CH]C(F)(F)F(12261)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {8,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
10 C u1 p0 c0 {6,S} {8,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {12,S}
12 H u0 p0 c0 {11,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1221.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.574334,0.111515,-0.000191587,1.70917e-07,-5.97137e-11,-146741,34.8521], Tmin=(100,'K'), Tmax=(802.216,'K')), NASAPolynomial(coeffs=[11.8578,0.034332,-1.88576e-05,3.76381e-09,-2.64897e-13,-148247,-19.3389], Tmin=(802.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1221.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H))"""),
)

species(
    label = 'FOC(F)(F)[CH][CH]C(F)(F)F(12262)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {8,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
10 C u1 p0 c0 {8,S} {11,S} {12,S}
11 C u1 p0 c0 {9,S} {10,S} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-885.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,253,522,585,610,849,1160,1215,1402,3000,3050,390,425,1340,1360,335,370,303.95,303.955,2126.24,2126.24],'cm^-1')),
        HinderedRotor(inertia=(0.386013,'amu*angstrom^2'), symmetry=1, barrier=(25.3068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115234,'amu*angstrom^2'), symmetry=1, barrier=(7.55481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115234,'amu*angstrom^2'), symmetry=1, barrier=(7.55484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386019,'amu*angstrom^2'), symmetry=1, barrier=(25.3068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0241268,0.0963569,-0.000132792,9.48721e-08,-2.72877e-11,-106322,35.6033], Tmin=(100,'K'), Tmax=(845.374,'K')), NASAPolynomial(coeffs=[13.7848,0.0310176,-1.68551e-05,3.44313e-09,-2.4945e-13,-108656,-28.697], Tmin=(845.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-885.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFF) + radical(CCJCO) + radical(Csj(Cs-CsHH)(Cs-F1sF1sF1s)(H))"""),
)

species(
    label = '[O]C(F)(F)[CH]C(F)C(F)(F)F(12251)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
10 C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1146.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,193,295,551,588,656,1146,1192,1350,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,335.095,335.246,335.364,2389.89],'cm^-1')),
        HinderedRotor(inertia=(0.0915993,'amu*angstrom^2'), symmetry=1, barrier=(7.30079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.091539,'amu*angstrom^2'), symmetry=1, barrier=(7.29547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205035,'amu*angstrom^2'), symmetry=1, barrier=(16.3554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.280406,0.101157,-0.000145734,1.07253e-07,-3.13452e-11,-137736,33.3409], Tmin=(100,'K'), Tmax=(837.562,'K')), NASAPolynomial(coeffs=[15.0411,0.0279849,-1.46882e-05,2.94523e-09,-2.10694e-13,-140303,-37.8605], Tmin=(837.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1146.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + radical(O2sj(Cs-CsF1sF1s)) + radical(CCJCO)"""),
)

species(
    label = 'FO[C](F)C(F)[CH]C(F)(F)F(12263)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {11,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {13,S}
11 C u1 p0 c0 {5,S} {7,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-853.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,395,473,707,1436,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.26294,'amu*angstrom^2'), symmetry=1, barrier=(6.0455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262954,'amu*angstrom^2'), symmetry=1, barrier=(6.04583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262994,'amu*angstrom^2'), symmetry=1, barrier=(6.04675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0460567,'amu*angstrom^2'), symmetry=1, barrier=(37.1087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.701497,0.11446,-0.000196213,1.74865e-07,-6.1224e-11,-102467,33.9671], Tmin=(100,'K'), Tmax=(793.757,'K')), NASAPolynomial(coeffs=[12.2067,0.0351937,-1.95519e-05,3.92268e-09,-2.7708e-13,-104068,-22.5045], Tmin=(793.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-853.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O][C](F)C(F)C(F)C(F)(F)F(12264)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  O u1 p2 c0 {11,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {7,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1137.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,259,529,569,1128,1321,1390,3140,193,295,551,588,656,1146,1192,1350,395,473,707,1436,180,180,180,2169.31],'cm^-1')),
        HinderedRotor(inertia=(0.332575,'amu*angstrom^2'), symmetry=1, barrier=(7.64656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333213,'amu*angstrom^2'), symmetry=1, barrier=(7.66121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332874,'amu*angstrom^2'), symmetry=1, barrier=(7.65344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.777206,0.11954,-0.000218587,1.99923e-07,-6.94082e-11,-136676,34.0636], Tmin=(100,'K'), Tmax=(856.637,'K')), NASAPolynomial(coeffs=[10.9164,0.0351099,-1.85191e-05,3.5761e-09,-2.44487e-13,-137585,-14.1541], Tmin=(856.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1137.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sCsH)(F1s)(O2s-H))"""),
)

species(
    label = 'FOC(F)(F)C(F)[CH][C](F)F(9381)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {7,S}
7  O u0 p2 c0 {6,S} {9,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u1 p0 c0 {8,S} {11,S} {13,S}
11 C u1 p0 c0 {4,S} {5,S} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-826.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,3025,407.5,1350,352.5,190,488,555,1236,1407,214.146,214.147,1658.75,1658.76],'cm^-1')),
        HinderedRotor(inertia=(0.961882,'amu*angstrom^2'), symmetry=1, barrier=(31.3011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280869,'amu*angstrom^2'), symmetry=1, barrier=(9.14005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280871,'amu*angstrom^2'), symmetry=1, barrier=(9.14002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961865,'amu*angstrom^2'), symmetry=1, barrier=(31.3011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.67483,0.113547,-0.0001924,1.7094e-07,-6.00529e-11,-99293.1,36.3622], Tmin=(100,'K'), Tmax=(779.764,'K')), NASAPolynomial(coeffs=[12.1722,0.0357837,-1.99928e-05,4.03098e-09,-2.85984e-13,-100936,-20.1091], Tmin=(779.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-826.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(F)(F)C(F)C(F)[C](F)F(5576)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
11 C u1 p0 c0 {5,S} {6,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1083.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,259,529,569,1128,1321,1390,3140,351,323,533,609,664,892,1120,1201,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.228149,'amu*angstrom^2'), symmetry=1, barrier=(5.2456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230334,'amu*angstrom^2'), symmetry=1, barrier=(5.29583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.90025,'amu*angstrom^2'), symmetry=1, barrier=(20.6985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750567,0.113602,-0.000187716,1.59767e-07,-5.35738e-11,-130172,34.4539], Tmin=(100,'K'), Tmax=(789.559,'K')), NASAPolynomial(coeffs=[14.0802,0.0309951,-1.65828e-05,3.28345e-09,-2.30469e-13,-132281,-32.1165], Tmin=(789.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1083.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    E0 = (-450.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-290.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-63.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-60.6101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-442.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-387.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-205.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-364.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-305.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-316.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-424.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-136.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-135.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-66.6833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-269.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-322.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-375.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-67.4763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-265.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-70.7872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-216.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-52.1855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-233.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    products = ['CF2O(49)', 'FC=CC(F)(F)F(1027)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    products = ['[O]C(F)(F)C([CH]F)C(F)(F)F(6991)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(F)(F)F(462)', '[O]C(F)(F)[CH]F(841)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'F[C](F)C(F)[CH]C(F)(F)F(1140)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    products = ['FC1C(C(F)(F)F)OC1(F)F(12107)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    products = ['OC(F)(F)C(F)=CC(F)(F)F(7536)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', '[O]C(F)(F)C=CC(F)(F)F(12255)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(65.2937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C](F)F(263)', 'FC=CC(F)(F)F(1027)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.88073e-05,'m^3/(mol*s)'), n=3.19085, Ea=(21.4794,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_N-4R!H->C_Sp-2R!H=1R!H_Ext-2R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_N-4R!H->C_Sp-2R!H=1R!H_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]C(F)(F)C(F)=CC(F)(F)F(12256)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.72227,'m^3/(mol*s)'), n=2.06948, Ea=(11.5395,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.15211669081533294, var=0.5516637502444228, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_N-Sp-7R!H=4CCClCl',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R_Ext-4CCl-R_N-Sp-7R!H=4CCClCl"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'O=C(F)C(F)[CH]C(F)(F)F(12257)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(33.3137,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CF2O(49)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(73.8291,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(F)(F)C(F)C=C(F)F(3661)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C](F)F(263)', 'F[CH][CH]C(F)(F)F(9334)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F2(78)', '[O]C(F)=C[CH]C(F)(F)F(12258)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(26.1815,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'O=C(F)[C](F)[CH]C(F)(F)F(12259)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(278.084,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    products = ['[O]C(F)(F)[C](F)CC(F)(F)F(12260)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    products = ['OC(F)(F)[C](F)[CH]C(F)(F)F(12261)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['FOC(F)(F)[CH][CH]C(F)(F)F(12262)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(116.497,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)(F)[CH]C(F)C(F)(F)F(12251)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(179.442,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FO[C](F)C(F)[CH]C(F)(F)F(12263)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(81.3093,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C](F)C(F)C(F)C(F)(F)F(12264)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(219.71,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FOC(F)(F)C(F)[CH][C](F)F(9381)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(73.5164,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(F)(F)C(F)C(F)[C](F)F(5576)'],
    products = ['[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(149.382,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3636',
    isomers = [
        '[O]C(F)(F)C(F)[CH]C(F)(F)F(12104)',
    ],
    reactants = [
        ('CF2O(49)', 'FC=CC(F)(F)F(1027)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3636',
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

