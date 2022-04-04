species(
    label = '[O]C(C(F)F)C1(F)[CH]O1(9867)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u1 p0 c0 {4,S} {6,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-483.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,2950,1000,180,180,386.539,1071.9,1071.96,1072.02],'cm^-1')),
        HinderedRotor(inertia=(0.559102,'amu*angstrom^2'), symmetry=1, barrier=(12.8548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00779046,'amu*angstrom^2'), symmetry=1, barrier=(40.3435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.859036,0.0718897,-8.01856e-05,4.64122e-08,-1.08614e-11,-58053,27.622], Tmin=(100,'K'), Tmax=(1028.75,'K')), NASAPolynomial(coeffs=[12.8288,0.0253485,-1.23245e-05,2.43564e-09,-1.74405e-13,-60515.8,-30.4644], Tmin=(1028.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-483.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(CC(C)OJ) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
)

species(
    label = 'FC1=CO1(3064)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-29.4903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,326,540,652,719,1357,788.854,788.882],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3010.76,'J/mol'), sigma=(4.86933,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=470.27 K, Pc=59.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.14414,0.0129866,6.60907e-06,-2.08037e-08,9.56931e-12,-3510.75,9.45209], Tmin=(100,'K'), Tmax=(957.696,'K')), NASAPolynomial(coeffs=[8.89006,0.00340757,-9.7292e-07,1.96224e-10,-1.667e-14,-5272.61,-21.4727], Tmin=(957.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.4903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + ring(Cyclopropene)"""),
)

species(
    label = 'O=CC(F)F(990)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-557.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([195,270,1147,1130,1359,1388,1409,3075,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.595284,'amu*angstrom^2'), symmetry=1, barrier=(13.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93683,0.00748989,0.000222302,-1.39306e-06,2.75418e-09,-67109.3,9.32225], Tmin=(10,'K'), Tmax=(169.939,'K')), NASAPolynomial(coeffs=[3.97443,0.0203448,-1.24424e-05,3.60772e-09,-4.02215e-13,-67130.4,8.62375], Tmin=(169.939,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-557.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODCC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]CC1(F)[CH]O1(9312)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {4,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-41.4931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,2850,1437.5,1250,1305,750,350,2950,1000,180,1536.43,1536.46,1536.58,1536.72],'cm^-1')),
        HinderedRotor(inertia=(0.198376,'amu*angstrom^2'), symmetry=1, barrier=(4.56105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3688.55,'J/mol'), sigma=(6.29887,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=576.14 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41455,0.0387057,-4.82418e-05,4.36624e-08,-1.70895e-11,-4937.09,19.4557], Tmin=(100,'K'), Tmax=(759.175,'K')), NASAPolynomial(coeffs=[3.55169,0.0268132,-1.30848e-05,2.55073e-09,-1.79548e-13,-4939.69,15.4029], Tmin=(759.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.4931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Cs(F)(C)-Cs-O2s) + radical(CCOJ) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
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
    label = '[O]C(F)C1(F)[CH]O1(9334)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
7 C u1 p0 c0 {3,S} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-258.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,391,562,707,872,1109,1210,1289,3137,2950,1000,180,1398.92,1399.24,1399.57],'cm^-1')),
        HinderedRotor(inertia=(0.301115,'amu*angstrom^2'), symmetry=1, barrier=(6.92323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3703.32,'J/mol'), sigma=(6.2169,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=578.45 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96847,0.0434771,-4.0118e-05,1.89811e-08,-3.63124e-12,-31076.1,21.4554], Tmin=(100,'K'), Tmax=(1244.61,'K')), NASAPolynomial(coeffs=[10.4473,0.0162269,-7.27581e-06,1.38919e-09,-9.75822e-14,-33186.7,-21.3055], Tmin=(1244.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(F)(C)-Cs-O2s) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
)

species(
    label = '[O]C(C(F)F)C1O[C]1F(9868)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-468.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,395,473,707,1436,180,180,180,180,804.158,804.577,805.141,805.753],'cm^-1')),
        HinderedRotor(inertia=(0.0781542,'amu*angstrom^2'), symmetry=1, barrier=(36.0782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00145887,'amu*angstrom^2'), symmetry=1, barrier=(4.19958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3798.02,'J/mol'), sigma=(6.08017,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.24 K, Pc=38.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958476,0.0679038,-6.89763e-05,3.57722e-08,-7.49929e-12,-56255,28.0756], Tmin=(100,'K'), Tmax=(1140.54,'K')), NASAPolynomial(coeffs=[13.3694,0.0243782,-1.17341e-05,2.31385e-09,-1.65547e-13,-59086.1,-33.4322], Tmin=(1140.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCsFFH) + ring(O2s-Cs(F)-Cs(C)) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
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
    label = 'FC(F)[CH]C1(F)[CH]O1(9276)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
7  C u1 p0 c0 {5,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-341.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,2950,1000,216.813,216.972,1216.07,1887.57,1888.04],'cm^-1')),
        HinderedRotor(inertia=(0.00358474,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910547,'amu*angstrom^2'), symmetry=1, barrier=(30.4603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44906,0.0552202,-5.08539e-05,2.38228e-08,-4.51938e-12,-40991.3,25.8248], Tmin=(100,'K'), Tmax=(1250.84,'K')), NASAPolynomial(coeffs=[12.188,0.0208789,-9.6723e-06,1.87429e-09,-1.3266e-13,-43677.9,-28.3882], Tmin=(1250.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(CCJCO) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
)

species(
    label = 'FC(F)C1OC2OC21F(9872)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-792.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742023,0.0632259,-5.07469e-05,1.84273e-08,-2.60234e-12,-95242.4,20.6086], Tmin=(100,'K'), Tmax=(1683.31,'K')), NASAPolynomial(coeffs=[20.6496,0.015921,-8.59411e-06,1.73316e-09,-1.23023e-13,-101945,-85.8011], Tmin=(1683.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-792.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'O=C(C(F)F)C1(F)CO1(9921)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
9  C u0 p0 c0 {5,D} {6,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-845.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30802,0.0636701,-6.48236e-05,3.52102e-08,-7.98107e-12,-101603,25.6349], Tmin=(100,'K'), Tmax=(1039.56,'K')), NASAPolynomial(coeffs=[10.3853,0.0287423,-1.44249e-05,2.88917e-09,-2.08209e-13,-103491,-18.5098], Tmin=(1039.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-845.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + ring(Cs-Cs(F)-O2s)"""),
)

species(
    label = '[O]C=C(F)C([O])C(F)F(9922)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-618.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,323,467,575,827,1418,3010,987.5,1337.5,450,1655,328.791,329.766,330.244,331.202],'cm^-1')),
        HinderedRotor(inertia=(0.257434,'amu*angstrom^2'), symmetry=1, barrier=(20.1322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260925,'amu*angstrom^2'), symmetry=1, barrier=(20.1338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366327,0.0771565,-8.40577e-05,4.42582e-08,-9.11507e-12,-74236.4,29.8977], Tmin=(100,'K'), Tmax=(1183.09,'K')), NASAPolynomial(coeffs=[17.9805,0.0176022,-8.54918e-06,1.70849e-09,-1.23644e-13,-78404.1,-58.0414], Tmin=(1183.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-618.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ)"""),
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
    label = '[O]C(C1=CO1)C(F)F(9923)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-261.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,2950,1000,180,180,180,456.324,675.492,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156203,'amu*angstrom^2'), symmetry=1, barrier=(3.59141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156203,'amu*angstrom^2'), symmetry=1, barrier=(3.59141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.599687,0.0651005,-5.87433e-05,1.87837e-08,2.94183e-14,-31321.2,26.8924], Tmin=(100,'K'), Tmax=(1007.8,'K')), NASAPolynomial(coeffs=[19.362,0.00933383,-3.5762e-06,7.03775e-10,-5.33196e-14,-36052.7,-68.4825], Tmin=(1007.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(CC(C)OJ)"""),
)

species(
    label = '[O][CH]C(F)F(1636)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-255.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.637092,'amu*angstrom^2'), symmetry=1, barrier=(14.648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0921,0.0467254,-7.94794e-05,7.13576e-08,-2.49602e-11,-30696.7,15.5487], Tmin=(100,'K'), Tmax=(816.59,'K')), NASAPolynomial(coeffs=[6.94864,0.0153515,-7.91621e-06,1.55898e-09,-1.09043e-13,-31237,-5.34886], Tmin=(816.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'F[C]1[CH]O1(9215)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u1 p0 c0 {2,S} {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (97.7552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,661.567,661.569,2935.1],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18921,0.0134531,-5.29464e-07,-1.06374e-08,5.3494e-12,11790.3,13.5601], Tmin=(100,'K'), Tmax=(993.343,'K')), NASAPolynomial(coeffs=[8.19763,0.00359113,-1.19992e-06,2.57126e-10,-2.11232e-14,10286.8,-13.1285], Tmin=(993.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.7552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'CHF2(82)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CC1(F)[CH]O1(9439)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 C u0 p0 c0 {3,D} {4,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-236.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,2782.5,750,1395,475,1775,1000,460.482,460.504,460.583],'cm^-1')),
        HinderedRotor(inertia=(0.179896,'amu*angstrom^2'), symmetry=1, barrier=(27.063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0451,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18077,0.033308,-1.89584e-05,-1.77368e-09,3.39657e-12,-28374.4,18.7218], Tmin=(100,'K'), Tmax=(1050.95,'K')), NASAPolynomial(coeffs=[11.7755,0.00993684,-4.36596e-06,8.73092e-10,-6.4647e-14,-31117.2,-31.4985], Tmin=(1050.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-236.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cs-Cs(F)-O2s) + radical(CCsJO)"""),
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
    label = 'O=C(C(F)F)C1(F)[CH]O1(9924)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {5,D} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-661.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,195,270,1147,1130,1359,1388,1409,3075,375,552.5,462.5,1710,2950,1000,506.869,506.873,506.874,506.874,506.875],'cm^-1')),
        HinderedRotor(inertia=(0.205175,'amu*angstrom^2'), symmetry=1, barrier=(18.8264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103264,'amu*angstrom^2'), symmetry=1, barrier=(18.8264,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02505,0.0682422,-7.77914e-05,4.46931e-08,-1.02965e-11,-79495,26.8503], Tmin=(100,'K'), Tmax=(1047.03,'K')), NASAPolynomial(coeffs=[13.2033,0.0217177,-1.11398e-05,2.25494e-09,-1.63612e-13,-82045.3,-32.4626], Tmin=(1047.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-661.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + ring(Cs-Cs(F)-O2s) + radical(CCsJO)"""),
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
    label = '[O]C(=C1[CH]O1)C(F)F(9925)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {7,D} {8,S}
7  C u0 p0 c0 {4,S} {5,S} {6,D}
8  C u1 p0 c0 {3,S} {6,S} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-344.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,2950,1000,536.682,536.683,536.683,536.683,536.683,536.683,536.683,536.684,536.689],'cm^-1')),
        HinderedRotor(inertia=(0.000585283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970167,0.05287,-2.48806e-05,-2.06118e-08,1.55871e-11,-41283.4,22.0674], Tmin=(100,'K'), Tmax=(944.221,'K')), NASAPolynomial(coeffs=[20.8322,0.002888,2.5509e-07,-4.41752e-11,-3.05946e-15,-46557,-80.6792], Tmin=(944.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(C=C(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = 'O=C([CH]F)C1(F)[CH]O1(9926)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,D} {5,S} {8,S}
7  C u1 p0 c0 {3,S} {5,S} {9,S}
8  C u1 p0 c0 {2,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-300.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,2950,1000,235,1215,1347,1486,3221,410.986,410.986,410.989,410.993,410.999],'cm^-1')),
        HinderedRotor(inertia=(0.000998024,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322081,'amu*angstrom^2'), symmetry=1, barrier=(38.6049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29365,0.0569421,-5.77263e-05,2.85545e-08,-5.5682e-12,-36087.9,24.2594], Tmin=(100,'K'), Tmax=(1240.53,'K')), NASAPolynomial(coeffs=[14.2374,0.0152064,-7.26189e-06,1.43513e-09,-1.02995e-13,-39299.4,-40.9769], Tmin=(1240.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + ring(Cs-Cs(F)-O2s) + radical(CCsJO) + radical(Csj(CO-CsO2d)(F1s)(H))"""),
)

species(
    label = 'O[C](C(F)F)C1(F)[CH]O1(9927)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-537.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.511846,0.084523,-0.000120086,8.4148e-08,-2.036e-11,-64508.1,28.5094], Tmin=(100,'K'), Tmax=(635.618,'K')), NASAPolynomial(coeffs=[11.5632,0.0271819,-1.35711e-05,2.64166e-09,-1.8487e-13,-66159.6,-21.7383], Tmin=(635.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(C2CsJOH) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
)

species(
    label = '[O][C](C(F)F)C1(F)CO1(9928)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-521.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.685407,0.0784467,-0.00010508,7.64136e-08,-2.25541e-11,-62620.2,27.3951], Tmin=(100,'K'), Tmax=(823.4,'K')), NASAPolynomial(coeffs=[10.9626,0.0285214,-1.41306e-05,2.777e-09,-1.96816e-13,-64312.7,-20.1896], Tmin=(823.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-521.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'OC([C](F)F)C1(F)[CH]O1(9929)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {4,S} {6,S} {11,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-512.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,316,385,515,654,689,1295,1380,1390,370,380,2900,435,2950,1000,190,488,555,1236,1407,354.606,354.609,354.611,2249.94,2249.94],'cm^-1')),
        HinderedRotor(inertia=(0.393314,'amu*angstrom^2'), symmetry=1, barrier=(35.0969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134061,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112063,'amu*angstrom^2'), symmetry=1, barrier=(9.99932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600121,0.0802496,-0.000109977,8.01424e-08,-2.34963e-11,-61473.1,29.7948], Tmin=(100,'K'), Tmax=(831.428,'K')), NASAPolynomial(coeffs=[11.7398,0.0266572,-1.32908e-05,2.61666e-09,-1.85547e-13,-63325.5,-21.8913], Tmin=(831.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = '[O]C([C](F)F)C1(F)CO1(9930)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-496.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940309,0.0717985,-8.40643e-05,5.26595e-08,-1.34873e-11,-59592.4,28.1047], Tmin=(100,'K'), Tmax=(939.348,'K')), NASAPolynomial(coeffs=[11.2891,0.0277293,-1.36908e-05,2.71337e-09,-1.94228e-13,-61536.6,-21.1746], Tmin=(939.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-496.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'FOC([C]1[CH]O1)C(F)F(9931)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
8  C u1 p0 c0 {4,S} {6,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-211.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,2950,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0552715,0.0907664,-0.000121033,7.01406e-08,-1.02647e-11,-25321.5,28.7685], Tmin=(100,'K'), Tmax=(702.753,'K')), NASAPolynomial(coeffs=[16.6898,0.0162033,-4.82483e-06,6.3788e-10,-3.20125e-14,-28156.3,-49.1505], Tmin=(702.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CCsJO)"""),
)

species(
    label = '[O]C([C]1OC1F)C(F)F(9913)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
9  C u1 p0 c0 {4,S} {6,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-511.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02886,0.0711122,-8.79759e-05,6.13646e-08,-1.78193e-11,-61432.4,28.1654], Tmin=(100,'K'), Tmax=(827.642,'K')), NASAPolynomial(coeffs=[9.28607,0.0312055,-1.56506e-05,3.10717e-09,-2.22025e-13,-62799.2,-10.1089], Tmin=(827.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-511.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCsFFH) + ring(O2s-Cs(F)-Cs(C)) + radical(CC(C)OJ) + radical(C2CsJO)"""),
)

species(
    label = 'F[CH]C(OF)C1(F)[CH]O1(9932)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {3,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {4,S} {6,S} {11,S}
9  C u1 p0 c0 {2,S} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-163.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,316,385,515,654,689,1295,1380,1390,370,380,2900,435,2950,1000,334,575,1197,1424,3202,187.688,188.888,189.189,189.378,190.774,2377.58],'cm^-1')),
        HinderedRotor(inertia=(0.0580053,'amu*angstrom^2'), symmetry=1, barrier=(35.8129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.99709e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.120318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40324,'amu*angstrom^2'), symmetry=1, barrier=(35.8126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263393,0.089196,-0.000138636,1.16433e-07,-3.92407e-11,-19523.8,29.6531], Tmin=(100,'K'), Tmax=(764.63,'K')), NASAPolynomial(coeffs=[10.9345,0.0298812,-1.54279e-05,3.03907e-09,-2.13774e-13,-21053.6,-18.298], Tmin=(764.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCsFHH) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
)

species(
    label = '[O]C([CH]F)C1(F)OC1F(9933)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {11,S}
9  C u1 p0 c0 {3,S} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-482.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,1380,1390,370,380,2900,435,261,493,600,1152,1365,1422,3097,334,575,1197,1424,3202,381.403,381.455,381.582,381.751],'cm^-1')),
        HinderedRotor(inertia=(0.00115814,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00115749,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.941321,0.071681,-8.0798e-05,4.82356e-08,-1.17955e-11,-57964.9,27.3765], Tmin=(100,'K'), Tmax=(979.774,'K')), NASAPolynomial(coeffs=[11.6252,0.0280625,-1.40186e-05,2.79624e-09,-2.01007e-13,-60058.4,-23.9486], Tmin=(979.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-482.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCsFHH) + ring(Cs(F)(C)-Cs-O2s) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-199.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (165.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (312.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-24.1283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (186.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-190.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-135.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-151.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (131.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-0.67936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-97.9105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-139.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-128.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (126.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-1.76828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-19.0143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-33.3603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-56.3567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-199.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-102.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (186.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-30.6254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (219.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (43.3597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['FC1=CO1(3064)', 'O=CC(F)F(990)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', '[O]CC1(F)[CH]O1(9312)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33604e-10,'m^3/(mol*s)'), n=4.72997, Ea=(126.614,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', '[O]C(F)C1(F)[CH]O1(9334)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.96165e-06,'m^3/(mol*s)'), n=3.30609, Ea=(148.557,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s',), comment="""Estimated from node CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(C(F)F)C1O[C]1F(9868)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'FC(F)[CH]C1(F)[CH]O1(9276)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['FC(F)C1OC2OC21F(9872)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['O=C(C(F)F)C1(F)CO1(9921)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C(F)C([O])C(F)F(9922)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[O]C(C1=CO1)C(F)F(9923)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(35.1005,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FC1=CO1(3064)', '[O][CH]C(F)F(1636)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[C]1[CH]O1(9215)', 'O=CC(F)F(990)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.66944e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(77.7499,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHF2(82)', 'O=CC1(F)[CH]O1(9439)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.43165,'m^3/(mol*s)'), n=1.58832, Ea=(69.5254,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.2743586997156029, var=2.1590031255329776, Tref=1000.0, N=359, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'O=C(C(F)F)C1(F)[CH]O1(9924)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(36.6168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]1[CH]O1(9215)', '[O][CH]C(F)F(1636)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', '[O]C(=C1[CH]O1)C(F)F(9925)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(339.035,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', 'O=C([CH]F)C1(F)[CH]O1(9926)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(278.408,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['O[C](C(F)F)C1(F)[CH]O1(9927)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['[O][C](C(F)F)C1(F)CO1(9928)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeO;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OC([C](F)F)C1(F)[CH]O1(9929)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.248e+06,'s^-1'), n=0.992, Ea=(28.4919,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 441 used for R3H_SS_Cs;C_rad_out_noH;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_noH;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 27.7 to 28.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['[O]C([C](F)F)C1(F)CO1(9930)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;Cs_H_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FOC([C]1[CH]O1)C(F)F(9931)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(113.801,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    products = ['[O]C([C]1OC1F)C(F)F(9913)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(168.405,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(OF)C1(F)[CH]O1(9932)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(98.7168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C([CH]F)C1(F)OC1F(9933)'],
    products = ['[O]C(C(F)F)C1(F)[CH]O1(9867)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(241.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #3298',
    isomers = [
        '[O]C(C(F)F)C1(F)[CH]O1(9867)',
    ],
    reactants = [
        ('FC1=CO1(3064)', 'O=CC(F)F(990)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3298',
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

