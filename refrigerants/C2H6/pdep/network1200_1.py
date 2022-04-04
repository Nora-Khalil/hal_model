species(
    label = '[CH2]OC([CH2])=CC(2407)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (125.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,354.463,354.528],'cm^-1')),
        HinderedRotor(inertia=(0.00134129,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134166,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164338,'amu*angstrom^2'), symmetry=1, barrier=(14.6524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164274,'amu*angstrom^2'), symmetry=1, barrier=(14.6525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15368,0.0536272,-2.70607e-05,-6.06453e-09,7.10178e-12,15162.3,24.0545], Tmin=(100,'K'), Tmax=(979.13,'K')), NASAPolynomial(coeffs=[14.0787,0.0207765,-7.29885e-06,1.29079e-09,-8.98209e-14,11674.9,-42.9122], Tmin=(979.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(C=COCJ)"""),
)

species(
    label = 'CH2O(19)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C=CC(1692)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.759585,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74634,0.0218191,8.22302e-06,-2.14762e-08,8.55596e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.61,'K')), NASAPolynomial(coeffs=[6.82084,0.0192337,-7.45616e-06,1.36535e-09,-9.53184e-14,16028,-10.4336], Tmin=(1025.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(S)(24)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56195e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.5454e-10,-1.62957e-14,50691.8,6.78382], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]OC([CH2])=C(1198)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {6,S} {7,S}
5  C u1 p0 c0 {1,S} {10,S} {11,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (161.173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,495.731,496.317],'cm^-1')),
        HinderedRotor(inertia=(0.124568,'amu*angstrom^2'), symmetry=1, barrier=(21.8929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124944,'amu*angstrom^2'), symmetry=1, barrier=(21.8915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00219389,'amu*angstrom^2'), symmetry=1, barrier=(21.9183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3077.81,'J/mol'), sigma=(5.44408,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.75 K, Pc=43.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84023,0.0375698,-4.75596e-06,-2.68503e-08,1.49076e-11,19471.6,19.0917], Tmin=(100,'K'), Tmax=(935.944,'K')), NASAPolynomial(coeffs=[13.8335,0.0113782,-2.94999e-06,4.76385e-10,-3.44637e-14,16128.7,-43.8397], Tmin=(935.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=COCJ)"""),
)

species(
    label = 'CH2(T)(17)',
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
    label = '[CH2]C([O])=CC(574)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {1,S} {3,D} {5,S}
5  C u1 p0 c0 {4,S} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (52.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,270.949],'cm^-1')),
        HinderedRotor(inertia=(0.126895,'amu*angstrom^2'), symmetry=1, barrier=(6.6168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126614,'amu*angstrom^2'), symmetry=1, barrier=(6.61684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96215,0.0395246,-2.15042e-05,-2.39262e-09,4.75075e-12,6363.74,18.5338], Tmin=(100,'K'), Tmax=(938.6,'K')), NASAPolynomial(coeffs=[10.4572,0.016349,-5.28591e-06,8.75308e-10,-5.83597e-14,4195.22,-24.9683], Tmin=(938.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]O[C]=CC(1205)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u1 p0 c0 {1,S} {10,S} {11,S}
5  C u1 p0 c0 {1,S} {3,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (247.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,3998.97],'cm^-1')),
        HinderedRotor(inertia=(0.842713,'amu*angstrom^2'), symmetry=1, barrier=(19.3756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843352,'amu*angstrom^2'), symmetry=1, barrier=(19.3903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0605996,'amu*angstrom^2'), symmetry=1, barrier=(19.4175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03825,0.0359608,-1.09154e-05,-1.24125e-08,7.65519e-12,29876.5,21.3147], Tmin=(100,'K'), Tmax=(982.224,'K')), NASAPolynomial(coeffs=[10.8658,0.0163224,-5.83388e-06,1.04525e-09,-7.33128e-14,27355.6,-25.1202], Tmin=(982.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COCJ) + radical(C=CJO)"""),
)

species(
    label = 'CC=C1CCO1(2408)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {2,S} {6,D}
6  C u0 p0 c0 {4,S} {5,D} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-91.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80743,0.0292899,5.46722e-05,-1.00717e-07,4.40636e-11,-10879.8,15.9789], Tmin=(100,'K'), Tmax=(910.558,'K')), NASAPolynomial(coeffs=[16.6138,0.0137771,-1.36656e-06,5.03551e-11,-4.56866e-15,-15629.5,-65.3413], Tmin=(910.558,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(2methyleneoxetane)"""),
)

species(
    label = '[CH2]O[C]1CC1C(2493)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (221.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12638,0.0555456,-4.14516e-05,1.61497e-08,-2.55055e-12,26788,22.261], Tmin=(100,'K'), Tmax=(1494.41,'K')), NASAPolynomial(coeffs=[13.3547,0.0228148,-8.59849e-06,1.49375e-09,-9.87539e-14,23133.1,-41.6461], Tmin=(1494.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(C2CsJOCs) + radical(CsJOCC2)"""),
)

species(
    label = '[CH2][C]1OCC1C(715)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {2,S} {6,S}
6  C u1 p0 c0 {5,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (225.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13423,0.0568868,-5.00846e-05,2.67587e-08,-5.73917e-12,27247.8,21.3121], Tmin=(100,'K'), Tmax=(1288.24,'K')), NASAPolynomial(coeffs=[9.43752,0.0253786,-6.72936e-06,8.71675e-10,-4.5833e-14,25583.6,-19.0052], Tmin=(1288.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1([CH]C)CO1(2494)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {4,S} {12,S}
6  C u1 p0 c0 {2,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (241.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261731,0.0670971,-6.51388e-05,3.39042e-08,-6.6421e-12,29185.8,24.672], Tmin=(100,'K'), Tmax=(1491.99,'K')), NASAPolynomial(coeffs=[13.9118,0.0175466,-2.29794e-06,5.28439e-12,1.3211e-14,26554.5,-41.8113], Tmin=(1491.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]=CC(1693)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u1 p0 c0 {4,S} {9,S} {10,S}
4  C u1 p0 c0 {2,D} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12014e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.76,'K')), NASAPolynomial(coeffs=[10.7464,0.0143241,-5.20138e-06,8.69083e-10,-5.48388e-14,40045.6,-31.3798], Tmin=(2065.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][O](117)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88411,-0.00363912,3.28554e-05,-4.13625e-08,1.59638e-11,23210.8,7.47975], Tmin=(100,'K'), Tmax=(933.052,'K')), NASAPolynomial(coeffs=[6.69328,0.000290109,8.61345e-07,-1.56334e-10,7.33637e-15,21991.3,-9.60391], Tmin=(933.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH2]OC(C)=[C]C(2495)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {1,S} {13,S} {14,S}
6  C u1 p0 c0 {3,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (204.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,357.024,357.034],'cm^-1')),
        HinderedRotor(inertia=(0.00132256,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106133,'amu*angstrom^2'), symmetry=1, barrier=(9.59991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106139,'amu*angstrom^2'), symmetry=1, barrier=(9.59995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106136,'amu*angstrom^2'), symmetry=1, barrier=(9.59991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22787,0.0567976,-4.52366e-05,1.90224e-08,-3.27222e-12,24647.5,23.5746], Tmin=(100,'K'), Tmax=(1367.6,'K')), NASAPolynomial(coeffs=[12.2501,0.0245593,-9.87723e-06,1.78565e-09,-1.21309e-13,21632.7,-33.052], Tmin=(1367.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=COCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=[C]C)OC(2496)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
3  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {13,S} {14,S}
6  C u1 p0 c0 {3,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
"""),
    E0 = (170.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,410.239,411.964],'cm^-1')),
        HinderedRotor(inertia=(0.0910233,'amu*angstrom^2'), symmetry=1, barrier=(11.0064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0924603,'amu*angstrom^2'), symmetry=1, barrier=(10.9842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904442,'amu*angstrom^2'), symmetry=1, barrier=(10.9975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918683,'amu*angstrom^2'), symmetry=1, barrier=(11.0061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896139,0.0615498,-5.4183e-05,2.52381e-08,-4.70671e-12,20681.8,22.9699], Tmin=(100,'K'), Tmax=(1293.99,'K')), NASAPolynomial(coeffs=[13.8449,0.0215221,-7.78221e-06,1.33204e-09,-8.79947e-14,17330.7,-42.8373], Tmin=(1293.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(C)O[CH2](1937)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u1 p0 c0 {4,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (117.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20711,0.0503686,-1.44483e-05,-1.98553e-08,1.18717e-11,14270.1,23.4543], Tmin=(100,'K'), Tmax=(985.136,'K')), NASAPolynomial(coeffs=[14.7358,0.0206139,-7.47743e-06,1.36924e-09,-9.78804e-14,10382.9,-47.8116], Tmin=(985.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]C=C([CH2])OC(711)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {6,S} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (84.6395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,989.931,990.243],'cm^-1')),
        HinderedRotor(inertia=(0.755733,'amu*angstrom^2'), symmetry=1, barrier=(17.3758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.078484,'amu*angstrom^2'), symmetry=1, barrier=(17.3763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75601,'amu*angstrom^2'), symmetry=1, barrier=(17.3822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755727,'amu*angstrom^2'), symmetry=1, barrier=(17.3756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05734,0.0529839,-1.59191e-05,-2.34164e-08,1.46719e-11,10296.4,22.1957], Tmin=(100,'K'), Tmax=(950.276,'K')), NASAPolynomial(coeffs=[16.2058,0.0178087,-5.52334e-06,9.50023e-10,-6.74823e-14,6126.51,-56.9063], Tmin=(950.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(38)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u2 p0 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82364,-0.000909745,3.21388e-05,-3.73491e-08,1.33095e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.803,'K')), NASAPolynomial(coeffs=[4.3048,0.0094308,-3.27565e-06,5.95137e-10,-4.2732e-14,40709.2,1.84239], Tmin=(960.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[C]=C(141)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 C u0 p0 c0 {4,D} {5,S} {8,S}
4 C u1 p0 c0 {1,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (283.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,607.789,608.037],'cm^-1')),
        HinderedRotor(inertia=(0.0541358,'amu*angstrom^2'), symmetry=1, barrier=(14.1943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0541755,'amu*angstrom^2'), symmetry=1, barrier=(14.1925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0632,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71888,0.0199752,1.113e-05,-3.28519e-08,1.53098e-11,34186,17.0661], Tmin=(100,'K'), Tmax=(932.137,'K')), NASAPolynomial(coeffs=[10.6431,0.00688566,-1.46263e-06,2.25522e-10,-1.75128e-14,31800.1,-25.4805], Tmin=(932.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OCC1C(533)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {2,S} {6,D}
6  C u0 p0 c0 {5,D} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-86.9896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79382,0.0259982,7.19633e-05,-1.2448e-07,5.40785e-11,-10361.9,14.8353], Tmin=(100,'K'), Tmax=(906.543,'K')), NASAPolynomial(coeffs=[18.8977,0.00966457,1.14233e-06,-4.42366e-10,2.86424e-14,-15892.9,-79.4051], Tmin=(906.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.9896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'C=CC(=C)OC(1931)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-101.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07036,0.0520215,-1.27095e-05,-2.70499e-08,1.5953e-11,-12065.3,19.3878], Tmin=(100,'K'), Tmax=(953.739,'K')), NASAPolynomial(coeffs=[16.5593,0.0171991,-5.34292e-06,9.34108e-10,-6.73975e-14,-16390.6,-61.7901], Tmin=(953.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[CH][C]1CCO1(2497)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {2,S} {6,S}
6  C u1 p0 c0 {4,S} {5,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (224.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02864,0.0369561,8.73484e-06,-3.86797e-08,1.96073e-11,27028.1,21.8341], Tmin=(100,'K'), Tmax=(851.751,'K')), NASAPolynomial(coeffs=[8.5683,0.0260749,-7.02518e-06,9.89799e-10,-5.91191e-14,25194.8,-12.8893], Tmin=(851.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'H(6)',
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
    label = '[CH2]OC(=C)C=C(1936)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  C u1 p0 c0 {1,S} {12,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (90.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3000,3100,440,815,1455,1000,360.403,360.45,360.459],'cm^-1')),
        HinderedRotor(inertia=(0.193067,'amu*angstrom^2'), symmetry=1, barrier=(17.8024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193148,'amu*angstrom^2'), symmetry=1, barrier=(17.8023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193077,'amu*angstrom^2'), symmetry=1, barrier=(17.8041,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1888,0.0498984,-1.42919e-05,-2.44687e-08,1.50517e-11,11023.2,21.4482], Tmin=(100,'K'), Tmax=(950.727,'K')), NASAPolynomial(coeffs=[16.5487,0.0144118,-4.27421e-06,7.42413e-10,-5.42788e-14,6785.72,-58.8038], Tmin=(950.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]CC(=C)O[CH2](2498)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {13,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (184.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,570.298,574.654,576.121],'cm^-1')),
        HinderedRotor(inertia=(0.00965765,'amu*angstrom^2'), symmetry=1, barrier=(2.21162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00969965,'amu*angstrom^2'), symmetry=1, barrier=(2.27231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.059562,'amu*angstrom^2'), symmetry=1, barrier=(14.0329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611043,'amu*angstrom^2'), symmetry=1, barrier=(14.0491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24537,0.0498588,-1.52172e-05,-1.83386e-08,1.12215e-11,22342,25.9293], Tmin=(100,'K'), Tmax=(986.575,'K')), NASAPolynomial(coeffs=[14.5014,0.0204182,-7.40871e-06,1.35573e-09,-9.67862e-14,18543.5,-43.8387], Tmin=(986.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=COCJ)"""),
)

species(
    label = '[CH]=C(CC)O[CH2](2499)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {1,S} {12,S} {13,S}
6  C u1 p0 c0 {4,D} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (226.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,453.953,454.359],'cm^-1')),
        HinderedRotor(inertia=(0.000817613,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0878056,'amu*angstrom^2'), symmetry=1, barrier=(12.8702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880767,'amu*angstrom^2'), symmetry=1, barrier=(12.8733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881042,'amu*angstrom^2'), symmetry=1, barrier=(12.872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13114,0.0526877,-2.157e-05,-1.30068e-08,9.65091e-12,27379.2,24.9558], Tmin=(100,'K'), Tmax=(984.877,'K')), NASAPolynomial(coeffs=[14.9695,0.0198428,-7.12194e-06,1.29458e-09,-9.21146e-14,23520.5,-47.3467], Tmin=(984.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=COCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=CC)OC(2500)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {1,S} {4,D} {6,S}
6  C u2 p0 c0 {5,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (144.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06485,0.0551733,-2.18386e-05,-8.80647e-09,6.76185e-12,17542.2,22.3088], Tmin=(100,'K'), Tmax=(1031.99,'K')), NASAPolynomial(coeffs=[12.6161,0.0291823,-1.13598e-05,2.05948e-09,-1.42579e-13,14157.9,-38.6287], Tmin=(1031.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-54.8131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (400.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (253.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (449.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-46.5288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (176.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (75.6316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (62.3565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (81.7611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (178.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (373.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (186.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (35.3297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (94.9741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (0.834066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (447.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-46.5288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-29.8399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (75.6316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (135.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (118.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (191.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (9.25568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['CH2O(19)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH2]OC([CH2])=C(1198)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.07568e+09,'m^3/(mol*s)'), n=-0.779019, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_4CbCdCt->Cd',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_4CbCdCt->Cd
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(17)', '[CH2]C([O])=CC(574)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(17)', '[CH2]O[C]=CC(1205)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['CC=C1CCO1(2408)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['[CH2]O[C]1CC1C(2493)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['[CH2][C]1OCC1C(715)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['[CH2]C1([CH]C)CO1(2494)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(117.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2O(19)', '[CH2][C]=CC(1693)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(19.7202,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][O](117)', 'C=C=CC(1692)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.01583e+10,'m^3/(mol*s)'), n=-1.23787, Ea=(20.0164,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08456738906112354, var=37.926887421795335, Tref=1000.0, N=24, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][O](117)', '[CH2][C]=CC(1693)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OC(C)=[C]C(2495)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=[C]C)OC(2496)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['[CH2]C=C(C)O[CH2](1937)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['[CH2]C=C([CH2])OC(711)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C(38)', '[CH2]O[C]=C(141)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['C=C1OCC1C(533)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['C=CC(=C)OC(1931)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OC([CH2])=CC(2407)'],
    products = ['C[CH][C]1CCO1(2497)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(6)', '[CH2]OC(=C)C=C(1936)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3407.27,'m^3/(mol*s)'), n=1.24181, Ea=(12.7419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.8400212858574567, var=0.41794239404157857, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_Sp-4CCl-2C_N-5C-inRing',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_N-4R!H->O_N-4BrCClFINPSSi-inRing_N-4BrCClFINPSSi->F_Ext-4CCl-R_N-Sp-5R!H#4CCCClClCl_5R!H->C_Sp-4CCl-2C_N-5C-inRing"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]CC(=C)O[CH2](2498)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(CC)O[CH2](2499)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]C(=CC)OC(2500)'],
    products = ['[CH2]OC([CH2])=CC(2407)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1200',
    isomers = [
        '[CH2]OC([CH2])=CC(2407)',
    ],
    reactants = [
        ('CH2O(19)', 'C=C=CC(1692)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1200',
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

