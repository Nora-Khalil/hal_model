species(
    label = '[O]C(F)(F)O[C](F)CF(10046)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-834.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,375,486,586,623,710,938,1161,1294,395,473,707,1436,198.895,198.9,198.902,199.03],'cm^-1')),
        HinderedRotor(inertia=(0.396541,'amu*angstrom^2'), symmetry=1, barrier=(11.1402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396691,'amu*angstrom^2'), symmetry=1, barrier=(11.1404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15742,'amu*angstrom^2'), symmetry=1, barrier=(32.5011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390834,0.0852281,-0.000125217,9.37337e-08,-2.77979e-11,-100209,27.5629], Tmin=(100,'K'), Tmax=(826.567,'K')), NASAPolynomial(coeffs=[13.2866,0.0228212,-1.19642e-05,2.38976e-09,-1.70323e-13,-102340,-32.1956], Tmin=(826.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-834.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs))"""),
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
    label = 'O=C(F)CF(476)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-611.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(15.7349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84069,0.0198378,-9.07231e-06,-3.50907e-10,9.31991e-13,-73601.2,9.53853], Tmin=(10,'K'), Tmax=(1222.46,'K')), NASAPolynomial(coeffs=[7.66923,0.0123956,-6.18005e-06,1.47455e-09,-1.37207e-13,-74917.2,-11.2551], Tmin=(1222.46,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-611.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C([O])(F)F(2271)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-372.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,464.782],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3646.55,'J/mol'), sigma=(6.01309,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.58 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92149,0.00463651,6.40183e-05,-1.39109e-07,8.57731e-11,-44748,10.2194], Tmin=(10,'K'), Tmax=(572.951,'K')), NASAPolynomial(coeffs=[4.84692,0.0160279,-1.25421e-05,4.35666e-09,-5.55356e-13,-45147.1,3.71309], Tmin=(572.951,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-372.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C([O])(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]CF(125)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p1 c0 {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-105.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,734,1109,1255,1358,2983,3011,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(0.0460479,'amu*angstrom^2'), symmetry=1, barrier=(35.5232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4418.31,'J/mol'), sigma=(4.687e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8659,0.0107238,1.79401e-05,-3.81644e-08,1.95617e-11,-12729.5,8.48308], Tmin=(10,'K'), Tmax=(672.698,'K')), NASAPolynomial(coeffs=[3.4045,0.0188633,-1.22417e-05,3.67097e-09,-4.17395e-13,-12789.5,9.6187], Tmin=(672.698,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-105.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]CF-2(467)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u2 p0 c0 {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-90.9927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([90,150,250,247,323,377,431,433,1065,1274,4000],'cm^-1')),
        HinderedRotor(inertia=(1.0536e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30451,0.015538,-7.22213e-06,9.68438e-10,5.66253e-14,-10919.1,11.642], Tmin=(100,'K'), Tmax=(1819.27,'K')), NASAPolynomial(coeffs=[8.39321,0.00776382,-3.62734e-06,6.82717e-10,-4.58694e-14,-13335.7,-17.5062], Tmin=(1819.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC[C](F)O[C](F)F(2217)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {4,S} {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-686.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,493,600,700,1144,1293,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.29083,'amu*angstrom^2'), symmetry=1, barrier=(6.68676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290858,'amu*angstrom^2'), symmetry=1, barrier=(6.68741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29081,'amu*angstrom^2'), symmetry=1, barrier=(6.6863,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641186,0.086932,-0.000168232,1.60804e-07,-5.74549e-11,-82452.4,26.987], Tmin=(100,'K'), Tmax=(860.066,'K')), NASAPolynomial(coeffs=[7.10834,0.0279478,-1.49451e-05,2.90845e-09,-1.9937e-13,-82495.7,2.97705], Tmin=(860.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-686.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs)) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'FCC1(F)OC(F)(F)O1(10327)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-1192.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56675,0.0585674,-5.53225e-05,2.70495e-08,-5.56067e-12,-143285,19.8748], Tmin=(100,'K'), Tmax=(1125.63,'K')), NASAPolynomial(coeffs=[10.1464,0.0280796,-1.46955e-05,2.98809e-09,-2.16765e-13,-145216,-22.5324], Tmin=(1125.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1192.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFOO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + ring(Cs-O2s-Cs(F)(F)-O2s)"""),
)

species(
    label = 'OC(F)(F)OC(F)=CF(5591)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1159.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,326,540,652,719,1357,194,682,905,1196,1383,3221,203.169,204.022,204.332],'cm^-1')),
        HinderedRotor(inertia=(0.567165,'amu*angstrom^2'), symmetry=1, barrier=(16.57,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356864,'amu*angstrom^2'), symmetry=1, barrier=(10.5178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934662,'amu*angstrom^2'), symmetry=1, barrier=(27.6895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.399677,0.0813805,-0.000104937,6.58641e-08,-1.61667e-11,-139382,25.0539], Tmin=(100,'K'), Tmax=(998.472,'K')), NASAPolynomial(coeffs=[16.3365,0.0175357,-9.02303e-06,1.82403e-09,-1.32258e-13,-142564,-51.8078], Tmin=(998.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1159.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = '[O][C](F)F(2186)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53788e-05,4.83907e-09,3.86776e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58999,0.00108969,-4.53602e-07,1.24872e-10,-1.13713e-14,-32130.4,-16.3888], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
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
    label = 'C=C(F)OC([O])(F)F(2323)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-726.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,326,540,652,719,1357,2950,3100,1380,975,1025,1650,182.812,326.948,326.955],'cm^-1')),
        HinderedRotor(inertia=(0.28278,'amu*angstrom^2'), symmetry=1, barrier=(21.5227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879175,'amu*angstrom^2'), symmetry=1, barrier=(21.5235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707644,0.0640665,-6.8219e-05,3.41441e-08,-6.54883e-12,-87260.2,24.0712], Tmin=(100,'K'), Tmax=(1286.12,'K')), NASAPolynomial(coeffs=[18.4847,0.00877763,-3.73591e-06,7.19019e-10,-5.1584e-14,-91832.9,-66.1664], Tmin=(1286.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-726.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(CdCFO) + group(Cds-CdsHH) + radical(O2sj(Cs-F1sF1sO2s))"""),
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
    label = '[O]C(F)(F)OC(F)=CF(10335)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-893.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180,669.996],'cm^-1')),
        HinderedRotor(inertia=(0.600778,'amu*angstrom^2'), symmetry=1, barrier=(13.8131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11272,'amu*angstrom^2'), symmetry=1, barrier=(25.5836,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600003,0.0733601,-8.84163e-05,5.05756e-08,-1.12057e-11,-107332,25.2273], Tmin=(100,'K'), Tmax=(1108.19,'K')), NASAPolynomial(coeffs=[17.3525,0.012892,-6.569e-06,1.33764e-09,-9.78957e-14,-111045,-57.3147], Tmin=(1108.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-893.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-F1sF1sO2s))"""),
)

species(
    label = 'O=C(F)O[C](F)CF(10195)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {4,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {5,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-835.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,482,664,788,1296,1923,329.403,329.41,329.422,329.465],'cm^-1')),
        HinderedRotor(inertia=(0.377157,'amu*angstrom^2'), symmetry=1, barrier=(29.0486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113528,'amu*angstrom^2'), symmetry=1, barrier=(8.74505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.505032,'amu*angstrom^2'), symmetry=1, barrier=(38.8948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46821,0.0608493,-7.31496e-05,4.67422e-08,-1.22984e-11,-100405,23.601], Tmin=(100,'K'), Tmax=(911.016,'K')), NASAPolynomial(coeffs=[9.80812,0.0242307,-1.28558e-05,2.61956e-09,-1.90196e-13,-101925,-15.8569], Tmin=(911.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-835.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COFOO) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O][C](F)CF(478)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-241.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.77787,'amu*angstrom^2'), symmetry=1, barrier=(17.8848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43603,0.0332746,-3.37965e-05,1.70828e-08,-3.40461e-12,-28931.2,16.3061], Tmin=(100,'K'), Tmax=(1215.45,'K')), NASAPolynomial(coeffs=[9.73799,0.00924406,-4.14003e-06,8.16317e-10,-5.88299e-14,-30706.2,-20.3463], Tmin=(1215.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)(F)OC(F)[CH]F(5114)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u1 p0 c0 {4,S} {7,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-830.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,375,486,586,623,710,938,1161,1294,334,575,1197,1424,3202,259.546,259.577,262.061,262.703],'cm^-1')),
        HinderedRotor(inertia=(0.443607,'amu*angstrom^2'), symmetry=1, barrier=(21.1627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127795,'amu*angstrom^2'), symmetry=1, barrier=(6.15854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810562,'amu*angstrom^2'), symmetry=1, barrier=(38.9664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.572284,0.0800059,-0.000104125,6.79e-08,-1.75929e-11,-99761.5,27.2906], Tmin=(100,'K'), Tmax=(940.183,'K')), NASAPolynomial(coeffs=[14.1831,0.0220987,-1.17381e-05,2.39007e-09,-1.73446e-13,-102321,-37.5344], Tmin=(940.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-830.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'OC(F)(F)O[C](F)[CH]F(10336)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {5,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-902.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335738,0.105178,-0.000187762,1.63641e-07,-5.454e-11,-108350,29.6091], Tmin=(100,'K'), Tmax=(850.093,'K')), NASAPolynomial(coeffs=[13.536,0.0228195,-1.22889e-05,2.38512e-09,-1.63385e-13,-110091,-31.4298], Tmin=(850.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-902.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[CH2][C](F)OC(F)(F)OF(10337)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {7,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {5,S} {9,S}
9  C u1 p0 c0 {8,S} {10,S} {11,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-579.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,435,565,619,662,854,1178,1396,395,473,707,1436,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.214996,'amu*angstrom^2'), symmetry=1, barrier=(4.94318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213143,'amu*angstrom^2'), symmetry=1, barrier=(4.90057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04607,'amu*angstrom^2'), symmetry=1, barrier=(24.0512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78653,'amu*angstrom^2'), symmetry=1, barrier=(41.0758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.117631,0.0998402,-0.0001735,1.52773e-07,-5.24406e-11,-69584.4,30.1725], Tmin=(100,'K'), Tmax=(806.424,'K')), NASAPolynomial(coeffs=[12.1712,0.0267322,-1.4909e-05,2.97859e-09,-2.09304e-13,-71171.2,-24.0201], Tmin=(806.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(CsCFHO) + group(Cs-CsHHH) + group(CsFFOO) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)OC([O])(F)F(10338)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {10,S} {11,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-873.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.672241,0.0779989,-0.000101317,6.64624e-08,-1.7382e-11,-104949,26.9227], Tmin=(100,'K'), Tmax=(930.359,'K')), NASAPolynomial(coeffs=[13.5627,0.0225772,-1.19614e-05,2.4329e-09,-1.76414e-13,-107347,-34.3359], Tmin=(930.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-873.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFFO) + group(Cs-CsHHH) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(Csj(Cs-F1sF1sO2s)(H)(H))"""),
)

species(
    label = 'FC[C](F)O[C](F)OF(10339)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u1 p0 c0 {2,S} {5,S} {7,S}
9  C u1 p0 c0 {3,S} {5,S} {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-527.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,482,586,761,1411,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00878807,0.0956058,-0.000160202,1.37053e-07,-4.5994e-11,-63357.1,29.7171], Tmin=(100,'K'), Tmax=(798.923,'K')), NASAPolynomial(coeffs=[12.5428,0.0251818,-1.35798e-05,2.68722e-09,-1.8823e-13,-65115.1,-26.4069], Tmin=(798.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFHOO) + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs)) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O][C](F)OC(F)(F)CF(10340)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-832.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.800459,0.0691673,-7.37465e-05,3.83209e-08,-7.86008e-12,-100024,28.2619], Tmin=(100,'K'), Tmax=(1180.37,'K')), NASAPolynomial(coeffs=[15.7996,0.018338,-9.15253e-06,1.83815e-09,-1.33025e-13,-103565,-46.5876], Tmin=(1180.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-832.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
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
    E0 = (-346.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (24.7384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (44.4135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-338.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-321.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-325.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-100.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-175.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-231.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-323.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-8.67169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (9.88217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-206.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-288.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (2.19684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-212.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (39.1334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-118.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(F)O[C](F)CF(10046)'],
    products = ['CF2O(49)', 'O=C(F)CF(476)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C([O])(F)F(2271)', 'F[C]CF-2(467)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', 'FC[C](F)O[C](F)F(2217)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)(F)O[C](F)CF(10046)'],
    products = ['FCC1(F)OC(F)(F)O1(10327)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(F)O[C](F)CF(10046)'],
    products = ['OC(F)(F)OC(F)=CF(5591)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C](F)F(2186)', 'O=C(F)CF(476)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(61.188,'m^3/(mol*s)'), n=1.29312, Ea=(54.1265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.28976280499384915, var=2.1569028208455543, Tref=1000.0, N=37, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'C=C(F)OC([O])(F)F(2323)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(65.1497,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', '[O]C(F)(F)OC(F)=CF(10335)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(50675.7,'m^3/(mol*s)'), n=0.856706, Ea=(18.0964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3149223058909152, var=1.98898808922157, Tref=1000.0, N=74, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'O=C(F)O[C](F)CF(10195)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(43.1142,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CF2O(49)', '[O][C](F)CF(478)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.30416e-16,'m^3/(mol*s)'), n=6.08142, Ea=(48.0038,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.42494241855262277, var=2.081425725945606, Tref=1000.0, N=380, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C](F)F(2186)', '[O][C](F)CF(478)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C([O])(F)F(2271)', 'F[C]CF(125)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(257654,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(F)OC(F)[CH]F(5114)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.30951e+09,'s^-1'), n=1.14834, Ea=(136.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(F)(F)O[C](F)CF(10046)'],
    products = ['OC(F)(F)O[C](F)[CH]F(10336)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.41845e+06,'s^-1'), n=1.36569, Ea=(57.7423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_1H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_1H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C](F)OC(F)(F)OF(10337)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(94.084,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)(F)O[C](F)CF(10046)'],
    products = ['[CH2]C(F)(F)OC([O])(F)F(10338)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(134.386,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['FC[C](F)O[C](F)OF(10339)'],
    products = ['[O]C(F)(F)O[C](F)CF(10046)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(79.218,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(F)(F)O[C](F)CF(10046)'],
    products = ['[O][C](F)OC(F)(F)CF(10340)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(227.418,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #3464',
    isomers = [
        '[O]C(F)(F)O[C](F)CF(10046)',
    ],
    reactants = [
        ('CF2O(49)', 'O=C(F)CF(476)'),
        ('[O]C([O])(F)F(2271)', 'F[C]CF(125)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3464',
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

