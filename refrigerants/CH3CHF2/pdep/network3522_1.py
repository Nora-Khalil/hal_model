species(
    label = '[CH2]C(F)C(F)C(=C)O(10093)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u0 p0 c0 {6,D} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-426.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,284,328,853,1146,1135,1297,3239,259,529,569,1128,1321,1390,3140,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,594.407],'cm^-1')),
        HinderedRotor(inertia=(0.829712,'amu*angstrom^2'), symmetry=1, barrier=(19.0767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.827802,'amu*angstrom^2'), symmetry=1, barrier=(19.0328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5765,'amu*angstrom^2'), symmetry=1, barrier=(13.2549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828084,'amu*angstrom^2'), symmetry=1, barrier=(19.0393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.150394,0.0843943,-9.00914e-05,4.93478e-08,-1.07968e-11,-51112.7,27.0478], Tmin=(100,'K'), Tmax=(1105.65,'K')), NASAPolynomial(coeffs=[16.0395,0.0269108,-1.2105e-05,2.32462e-09,-1.64273e-13,-54626.3,-51.2036], Tmin=(1105.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-426.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C[CH]F(1341)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u1 p0 c0 {1,S} {2,S} {8,S}
4 C u0 p0 c0 {2,D} {6,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-39.4934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.618875,'amu*angstrom^2'), symmetry=1, barrier=(28.5117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2734.35,'J/mol'), sigma=(4.78057,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=427.10 K, Pc=56.79 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95088,0.00313686,8.55963e-05,-1.7635e-07,1.14205e-10,-4747,9.3909], Tmin=(10,'K'), Tmax=(483.604,'K')), NASAPolynomial(coeffs=[2.28776,0.0264346,-1.62621e-05,4.86442e-09,-5.65127e-13,-4697.72,15.0527], Tmin=(483.604,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-39.4934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CDC[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2F-CH2(62)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-75.6489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0750522,'amu*angstrom^2'), symmetry=1, barrier=(7.14957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0516,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2595.78,'J/mol'), sigma=(4.583,'angstroms'), dipoleMoment=(2,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.99845,0.00714341,1.60323e-05,-2.0007e-08,6.54682e-12,-9101.39,8.71379], Tmin=(10,'K'), Tmax=(1078.03,'K')), NASAPolynomial(coeffs=[3.84792,0.0159784,-7.77703e-06,1.83859e-09,-1.70814e-13,-9549.86,7.22073], Tmin=(1078.03,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-75.6489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[CH2]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(O)[C]F(661)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 C u0 p1 c0 {1,S} {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-45.5134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,315,622,1128],'cm^-1')),
        HinderedRotor(inertia=(0.949608,'amu*angstrom^2'), symmetry=1, barrier=(21.8334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95311,'amu*angstrom^2'), symmetry=1, barrier=(21.9139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04374,0.0341451,-1.44819e-05,-1.42936e-08,1.01491e-11,-5395.46,15.8839], Tmin=(100,'K'), Tmax=(949.392,'K')), NASAPolynomial(coeffs=[14.4867,0.00352083,-5.41025e-07,1.03936e-10,-1.11499e-14,-8740.61,-48.6744], Tmin=(949.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.5134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(CJ2_singlet-FC)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = '[CH2]C(F)C(=C)O(673)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {12,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u1 p0 c0 {3,S} {8,S} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
"""),
    E0 = (-204.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.796472,'amu*angstrom^2'), symmetry=1, barrier=(18.3125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796471,'amu*angstrom^2'), symmetry=1, barrier=(18.3124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796725,'amu*angstrom^2'), symmetry=1, barrier=(18.3183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0882,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988331,0.0602142,-5.81068e-05,2.80545e-08,-5.3222e-12,-24497.3,21.8994], Tmin=(100,'K'), Tmax=(1282.39,'K')), NASAPolynomial(coeffs=[15.1645,0.015996,-6.3852e-06,1.16632e-09,-8.03958e-14,-28133.2,-50.0186], Tmin=(1282.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCdH)(H)(H))"""),
)

species(
    label = 'C=C(O)C(F)C[CH]F(10326)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {2,S} {5,S} {12,S}
8  C u0 p0 c0 {6,D} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-431.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.260243,0.0848863,-9.51691e-05,5.60667e-08,-1.33277e-11,-51822.6,27.0655], Tmin=(100,'K'), Tmax=(1015.96,'K')), NASAPolynomial(coeffs=[14.2565,0.0297802,-1.38074e-05,2.67692e-09,-1.89768e-13,-54666.4,-40.6796], Tmin=(1015.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-431.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(F)C=C(O)CF(10327)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {12,S}
8  C u1 p0 c0 {4,S} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-416.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.68065,'amu*angstrom^2'), symmetry=1, barrier=(15.6495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682252,'amu*angstrom^2'), symmetry=1, barrier=(15.6863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0642607,'amu*angstrom^2'), symmetry=1, barrier=(15.6589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680023,'amu*angstrom^2'), symmetry=1, barrier=(15.6351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.427326,0.0809716,-8.53512e-05,4.71372e-08,-1.05758e-11,-49924.7,27.3742], Tmin=(100,'K'), Tmax=(1068.28,'K')), NASAPolynomial(coeffs=[13.9433,0.0303648,-1.4295e-05,2.79533e-09,-1.99171e-13,-52812.6,-38.7256], Tmin=(1068.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-416.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + group(Cs-CsHHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Csj(Cs-F1sCdH)(H)(H))"""),
)

species(
    label = '[CH2]C(F)CC(O)=CF(10328)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {15,S}
4  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {6,D} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-425.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410044,0.0811535,-8.70509e-05,4.94219e-08,-1.13872e-11,-51037.6,28.039], Tmin=(100,'K'), Tmax=(1043.51,'K')), NASAPolynomial(coeffs=[13.6876,0.0302581,-1.38917e-05,2.68317e-09,-1.89864e-13,-53808.7,-36.5833], Tmin=(1043.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(O)C(F)[CH]F(678)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {12,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {2,S} {4,S} {9,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-391.649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.897151,'amu*angstrom^2'), symmetry=1, barrier=(20.6273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.897634,'amu*angstrom^2'), symmetry=1, barrier=(20.6384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895062,'amu*angstrom^2'), symmetry=1, barrier=(20.5792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.642134,0.0718123,-8.22225e-05,4.69737e-08,-1.04976e-11,-46981.7,23.5994], Tmin=(100,'K'), Tmax=(1096.52,'K')), NASAPolynomial(coeffs=[15.6811,0.0169505,-7.17213e-06,1.34335e-09,-9.39469e-14,-50279.8,-50.3405], Tmin=(1096.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-391.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCdH)(F1s)(H))"""),
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
    label = '[CH]C(F)C(F)C(=C)O(10329)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {14,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  C u2 p0 c0 {5,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-184.646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,284,328,853,1146,1135,1297,3239,250,417,511,1155,1315,1456,3119,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105085,0.0842385,-9.34194e-05,5.19531e-08,-1.13956e-11,-22066.3,26.8083], Tmin=(100,'K'), Tmax=(1110.56,'K')), NASAPolynomial(coeffs=[17.2295,0.0225593,-1.01104e-05,1.94233e-09,-1.37508e-13,-25869.8,-57.6027], Tmin=(1110.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'O[C]1CCC(F)C1F(10330)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {8,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {12,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-487.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924398,0.0571042,-2.36742e-05,-9.73875e-09,7.62406e-12,-58457.2,24.8416], Tmin=(100,'K'), Tmax=(1029.76,'K')), NASAPolynomial(coeffs=[14.3877,0.026353,-1.02648e-05,1.89827e-09,-1.33883e-13,-62372.3,-46.0527], Tmin=(1029.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-487.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + ring(Cyclopentane) + radical(C2CsJOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C1(O)CC(F)C1F(10331)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {15,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
8  C u1 p0 c0 {4,S} {13,S} {14,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310418,0.0726764,-6.31623e-05,2.78844e-08,-4.89597e-12,-47366.9,25.7104], Tmin=(100,'K'), Tmax=(1370.12,'K')), NASAPolynomial(coeffs=[17.1614,0.0234807,-9.30313e-06,1.6778e-09,-1.14172e-13,-51984.4,-60.8923], Tmin=(1370.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs-Cs) + radical(CJC(C)2O) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)C(F)C(C)=O(10332)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,D} {4,S} {6,S}
8  C u1 p0 c0 {5,S} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-458.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.875741,0.0747488,-7.93962e-05,4.9148e-08,-1.30502e-11,-55047.1,26.6181], Tmin=(100,'K'), Tmax=(888.499,'K')), NASAPolynomial(coeffs=[8.84351,0.0388778,-1.88367e-05,3.70801e-09,-2.6446e-13,-56463,-10.8797], Tmin=(888.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-458.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = 'C=CC(F)C(=C)O(10333)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-308.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05757,'amu*angstrom^2'), symmetry=1, barrier=(24.3155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05783,'amu*angstrom^2'), symmetry=1, barrier=(24.3216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05512,'amu*angstrom^2'), symmetry=1, barrier=(24.2593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.672971,0.0675674,-6.127e-05,2.82164e-08,-5.1773e-12,-36991.9,23.3034], Tmin=(100,'K'), Tmax=(1310.61,'K')), NASAPolynomial(coeffs=[15.6041,0.021997,-9.11399e-06,1.68602e-09,-1.16573e-13,-40905.6,-52.7693], Tmin=(1310.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH2CHF(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-153.05,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0219,'amu')),
        NonlinearRotor(inertia=([7.59478,47.6085,55.2033],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([482.184,740.114,880.476,949.569,983.363,1189.2,1343.65,1421.6,1725.69,3171.53,3191.61,3269.97],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09164,-0.0073724,7.45741e-05,-1.12982e-07,5.61696e-11,-18407.3,6.78145], Tmin=(10,'K'), Tmax=(619.705,'K')), NASAPolynomial(coeffs=[1.44203,0.0189088,-1.12569e-05,3.25441e-09,-3.64262e-13,-18255.1,16.8744], Tmin=(619.705,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-153.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(O)[CH]F(676)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {9,S}
3 C u0 p0 c0 {2,S} {4,S} {5,D}
4 C u1 p0 c0 {1,S} {3,S} {6,S}
5 C u0 p0 c0 {3,D} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-212.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.44574,'amu*angstrom^2'), symmetry=1, barrier=(33.2404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44148,'amu*angstrom^2'), symmetry=1, barrier=(33.1425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0616,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85729,0.00846381,0.000121446,-2.61705e-07,1.62019e-10,-25556.5,9.95283], Tmin=(10,'K'), Tmax=(567.87,'K')), NASAPolynomial(coeffs=[5.43441,0.0299714,-2.15197e-05,7.2771e-09,-9.25586e-13,-26261.5,-1.39388], Tmin=(567.87,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-212.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CDC(O)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C(O)C(F)C(=C)F(9278)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {14,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-506.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,323,467,575,827,1418,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02862,'amu*angstrom^2'), symmetry=1, barrier=(23.6499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03093,'amu*angstrom^2'), symmetry=1, barrier=(23.703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03022,'amu*angstrom^2'), symmetry=1, barrier=(23.6869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386946,0.0777574,-8.20007e-05,4.42042e-08,-9.47243e-12,-60819.6,24.7716], Tmin=(100,'K'), Tmax=(1131.75,'K')), NASAPolynomial(coeffs=[15.7031,0.0236237,-1.02516e-05,1.9391e-09,-1.36023e-13,-64286.4,-51.0154], Tmin=(1131.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCsCdF) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(O)=CC([CH2])F(10334)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {2,S} {4,D} {7,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-101.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.628792,'amu*angstrom^2'), symmetry=1, barrier=(14.4572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629005,'amu*angstrom^2'), symmetry=1, barrier=(14.4621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.03678,'amu*angstrom^2'), symmetry=1, barrier=(14.4737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339231,'amu*angstrom^2'), symmetry=1, barrier=(14.4736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311952,0.0754171,-7.88495e-05,4.16608e-08,-8.6186e-12,-12067.1,27.3313], Tmin=(100,'K'), Tmax=(1183.47,'K')), NASAPolynomial(coeffs=[17.0504,0.0188436,-7.1458e-06,1.26952e-09,-8.63271e-14,-16029,-56.2418], Tmin=(1183.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Csj(Cs-F1sCdH)(H)(H)) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][CH]C(F)C(=C)O(10335)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u1 p0 c0 {3,S} {7,S} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-35.4614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.000507363,'amu*angstrom^2'), symmetry=1, barrier=(5.76063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973415,'amu*angstrom^2'), symmetry=1, barrier=(22.3807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.972202,'amu*angstrom^2'), symmetry=1, barrier=(22.3528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0427845,'amu*angstrom^2'), symmetry=1, barrier=(22.4181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99928,0.0671259,-6.43502e-05,3.26393e-08,-6.79315e-12,-4157.92,26.6701], Tmin=(100,'K'), Tmax=(1139.78,'K')), NASAPolynomial(coeffs=[12.1705,0.0279215,-1.27561e-05,2.46184e-09,-1.74083e-13,-6704.49,-28.6863], Tmin=(1139.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.4614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'OH(7)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(F)C(F)[C]=C(7449)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {7,D} {12,S} {13,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (4.69768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.268517,'amu*angstrom^2'), symmetry=1, barrier=(6.17374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268393,'amu*angstrom^2'), symmetry=1, barrier=(6.17088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12765,'amu*angstrom^2'), symmetry=1, barrier=(25.9268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3088.83,'J/mol'), sigma=(5.37861,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=482.47 K, Pc=45.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00568,0.0715471,-9.30537e-05,7.18406e-08,-2.32715e-11,667.563,25.3397], Tmin=(100,'K'), Tmax=(744.481,'K')), NASAPolynomial(coeffs=[8.08538,0.0335097,-1.64166e-05,3.2154e-09,-2.27404e-13,-386.601,-6.7269], Tmin=(744.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.69768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCsH)(H)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-HH))"""),
)

species(
    label = '[CH2]C(F)C(F)C(=C)[O](10336)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u0 p0 c0 {6,D} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-288.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,259,529,569,1128,1321,1390,3140,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,336.123,336.129,336.134],'cm^-1')),
        HinderedRotor(inertia=(0.0872712,'amu*angstrom^2'), symmetry=1, barrier=(6.99683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0872713,'amu*angstrom^2'), symmetry=1, barrier=(6.99686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307523,'amu*angstrom^2'), symmetry=1, barrier=(24.6558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607212,0.0792052,-9.27504e-05,5.94547e-08,-1.56226e-11,-34559.3,26.9201], Tmin=(100,'K'), Tmax=(916.895,'K')), NASAPolynomial(coeffs=[11.5015,0.031677,-1.49939e-05,2.91686e-09,-2.06576e-13,-36557,-24.6929], Tmin=(916.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2][CH]F(233)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (116.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0021411,'amu*angstrom^2'), symmetry=1, barrier=(6.24533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56892,0.0121147,-4.46308e-06,3.53003e-10,6.36796e-14,13988.6,11.2355], Tmin=(100,'K'), Tmax=(2106.76,'K')), NASAPolynomial(coeffs=[7.95462,0.00674548,-2.74609e-06,4.76054e-10,-2.99987e-14,11484.3,-14.7487], Tmin=(2106.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(Cs_P)"""),
)

species(
    label = 'C=[C]O(137)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {6,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16242,0.0134243,5.5644e-06,-1.95524e-08,9.36428e-12,12455.2,10.1543], Tmin=(100,'K'), Tmax=(925.608,'K')), NASAPolynomial(coeffs=[8.19869,0.00453472,-8.93509e-07,1.26097e-10,-9.46636e-15,10971.4,-16.7327], Tmin=(925.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(F)[CH]F(252)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 C u1 p0 c0 {2,S} {3,S} {9,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-106.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,3000,3100,440,815,1455,1000,334,575,1197,1424,3202,1191.69],'cm^-1')),
        HinderedRotor(inertia=(0.266097,'amu*angstrom^2'), symmetry=1, barrier=(6.1181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265589,'amu*angstrom^2'), symmetry=1, barrier=(6.10641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.0606,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18615,0.0452068,-6.92125e-05,6.54885e-08,-2.48594e-11,-12751.2,18.7364], Tmin=(100,'K'), Tmax=(793.498,'K')), NASAPolynomial(coeffs=[4.32719,0.0247775,-1.23772e-05,2.433e-09,-1.71221e-13,-12787.6,10.8139], Tmin=(793.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(H)(H)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(F)[C](F)C(=C)O(10337)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {14,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
5  C u1 p0 c0 {2,S} {4,S} {6,S}
6  C u0 p0 c0 {3,S} {5,S} {8,D}
7  C u1 p0 c0 {4,S} {10,S} {11,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-288.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,346,659,817,1284,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.970844,'amu*angstrom^2'), symmetry=1, barrier=(22.3216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970826,'amu*angstrom^2'), symmetry=1, barrier=(22.3212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970815,'amu*angstrom^2'), symmetry=1, barrier=(22.3209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970947,'amu*angstrom^2'), symmetry=1, barrier=(22.324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0503988,0.0840677,-9.1763e-05,5.00056e-08,-1.07139e-11,-34515.8,27.5841], Tmin=(100,'K'), Tmax=(1138.39,'K')), NASAPolynomial(coeffs=[17.8306,0.0215923,-9.4417e-06,1.79598e-09,-1.26585e-13,-38563.9,-60.4995], Tmin=(1138.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCdCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2][C](F)C(F)C(=C)O(10338)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {14,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {2,S} {4,S} {8,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  C u1 p0 c0 {6,S} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
"""),
    E0 = (-235.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,350,440,435,1725,212,367,445,1450,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,351.707],'cm^-1')),
        HinderedRotor(inertia=(0.110062,'amu*angstrom^2'), symmetry=1, barrier=(9.66111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167573,'amu*angstrom^2'), symmetry=1, barrier=(14.7102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110051,'amu*angstrom^2'), symmetry=1, barrier=(9.66115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301974,'amu*angstrom^2'), symmetry=1, barrier=(26.5358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269348,0.0874107,-0.000114161,7.98349e-08,-2.25159e-11,-28213.3,28.6421], Tmin=(100,'K'), Tmax=(862.666,'K')), NASAPolynomial(coeffs=[12.6283,0.0301048,-1.45179e-05,2.83068e-09,-2.00081e-13,-30345.6,-29.157], Tmin=(862.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-235.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH]=C(O)C(F)C([CH2])F(10339)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {13,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u1 p0 c0 {6,D} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-179.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,284,328,853,1146,1135,1297,3239,259,529,569,1128,1321,1390,3140,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.740641,'amu*angstrom^2'), symmetry=1, barrier=(17.0288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.89418,'amu*angstrom^2'), symmetry=1, barrier=(20.559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740644,'amu*angstrom^2'), symmetry=1, barrier=(17.0289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740721,'amu*angstrom^2'), symmetry=1, barrier=(17.0306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.107606,0.0874902,-0.000105013,6.42264e-08,-1.55711e-11,-21394.2,27.6092], Tmin=(100,'K'), Tmax=(1005.63,'K')), NASAPolynomial(coeffs=[15.8714,0.0247885,-1.14879e-05,2.22564e-09,-1.57765e-13,-24564.7,-48.5307], Tmin=(1005.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCsH)(H)(H)) + radical(Cds_P)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C=CC(=C)O(10340)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {5,D}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {6,S} {7,S}
5  C u0 p0 c0 {2,D} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-37.9415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.20546,'amu*angstrom^2'), symmetry=1, barrier=(27.7159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19709,'amu*angstrom^2'), symmetry=1, barrier=(27.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20216,'amu*angstrom^2'), symmetry=1, barrier=(27.64,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24334,0.0466737,7.91805e-07,-4.47883e-08,2.38854e-11,-4451.11,20.4972], Tmin=(100,'K'), Tmax=(920.207,'K')), NASAPolynomial(coeffs=[17.5625,0.0120976,-2.1178e-06,2.59876e-10,-1.90781e-14,-8994.01,-65.2413], Tmin=(920.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.9415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=C(F)C(=C)O(10341)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {6,D}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {7,S} {8,S}
6  C u0 p0 c0 {3,D} {9,S} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-233.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,280,518,736,852,873,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.995401,'amu*angstrom^2'), symmetry=1, barrier=(22.8862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993747,'amu*angstrom^2'), symmetry=1, barrier=(22.8482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86966,'amu*angstrom^2'), symmetry=1, barrier=(42.9872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811038,0.0596432,-3.50773e-05,-4.90938e-09,8.35882e-12,-27903.1,22.618], Tmin=(100,'K'), Tmax=(960.686,'K')), NASAPolynomial(coeffs=[17.1751,0.0160107,-5.2079e-06,9.11634e-10,-6.48198e-14,-32177.9,-61.5576], Tmin=(960.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(F)=CC(=C)O(10342)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {6,D}
4  C u0 p0 c0 {3,S} {5,D} {8,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u0 p0 c0 {3,D} {9,S} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-246.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,271,519,563,612,1379,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.09639,'amu*angstrom^2'), symmetry=1, barrier=(25.2081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09663,'amu*angstrom^2'), symmetry=1, barrier=(25.2137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09655,'amu*angstrom^2'), symmetry=1, barrier=(25.2119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709651,0.0613156,-3.55553e-05,-9.41959e-09,1.16373e-11,-29511.2,22.51], Tmin=(100,'K'), Tmax=(925.573,'K')), NASAPolynomial(coeffs=[18.5929,0.0126026,-2.9143e-06,4.21589e-10,-2.91786e-14,-34045.5,-68.9943], Tmin=(925.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C(O)C(F)[C](C)F(10343)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-442.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,174,267,591,721,1107,1278,1348,3273,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,350,440,435,1725,2950,3100,1380,975,1025,1650,296.092],'cm^-1')),
        HinderedRotor(inertia=(0.224838,'amu*angstrom^2'), symmetry=1, barrier=(13.9878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224837,'amu*angstrom^2'), symmetry=1, barrier=(13.9878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224838,'amu*angstrom^2'), symmetry=1, barrier=(13.9878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224838,'amu*angstrom^2'), symmetry=1, barrier=(13.9878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416916,0.082562,-9.23153e-05,5.54471e-08,-1.35887e-11,-53081,26.7052], Tmin=(100,'K'), Tmax=(981.461,'K')), NASAPolynomial(coeffs=[12.8337,0.0319563,-1.49725e-05,2.91103e-09,-2.06478e-13,-55518.3,-32.9663], Tmin=(981.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-442.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'C=C(O)[C](F)C(C)F(10344)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {4,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-494.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0850906,0.0825017,-8.11152e-05,3.97363e-08,-7.64771e-12,-59371.2,26.6658], Tmin=(100,'K'), Tmax=(1264.05,'K')), NASAPolynomial(coeffs=[19.1489,0.0216375,-8.89058e-06,1.645e-09,-1.14187e-13,-64233.8,-70.6344], Tmin=(1264.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-494.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsCdCsF1s)"""),
)

species(
    label = 'C=C([O])C(F)C(C)F(10345)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {14,S} {15,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-495.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.699042,0.0750912,-7.38031e-05,3.92356e-08,-8.63253e-12,-59424.7,25.1775], Tmin=(100,'K'), Tmax=(1077.7,'K')), NASAPolynomial(coeffs=[12.0667,0.032899,-1.5078e-05,2.90834e-09,-2.05517e-13,-61874.9,-30.5156], Tmin=(1077.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-495.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(F)C(C)F(10346)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-385.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,250,417,511,1155,1315,1456,3119,284,328,853,1146,1135,1297,3239,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.769465,'amu*angstrom^2'), symmetry=1, barrier=(17.6915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.770065,'amu*angstrom^2'), symmetry=1, barrier=(17.7053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769657,'amu*angstrom^2'), symmetry=1, barrier=(17.6959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769715,'amu*angstrom^2'), symmetry=1, barrier=(17.6973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101729,0.0844937,-8.9811e-05,4.86526e-08,-1.04834e-11,-46255.4,26.2193], Tmin=(100,'K'), Tmax=(1124.74,'K')), NASAPolynomial(coeffs=[16.6321,0.0257052,-1.14078e-05,2.18045e-09,-1.53802e-13,-49973.8,-55.4734], Tmin=(1124.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O)=CC(F)CF(10347)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {3,S} {6,D} {8,S}
8  C u1 p0 c0 {7,S} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-482.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0484634,0.0825437,-8.4803e-05,4.38371e-08,-8.88307e-12,-57878,28.1551], Tmin=(100,'K'), Tmax=(1205.94,'K')), NASAPolynomial(coeffs=[18.5256,0.0209357,-8.17287e-06,1.47483e-09,-1.01138e-13,-62357.9,-64.9323], Tmin=(1205.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-482.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C(O)C(F)[CH]CF(10348)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u1 p0 c0 {4,S} {5,S} {12,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
"""),
    E0 = (-429.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (121.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.585327,0.078242,-8.04648e-05,4.38783e-08,-9.81907e-12,-51541.7,27.3137], Tmin=(100,'K'), Tmax=(1064.5,'K')), NASAPolynomial(coeffs=[12.9212,0.031888,-1.51468e-05,2.97139e-09,-2.11976e-13,-54168,-32.971], Tmin=(1064.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-429.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S)"""),
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
    E0 = (-66.7633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (95.2955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (269.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-111.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (43.4394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-109.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (145.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (183.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-213.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-207.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-65.3069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-17.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-178.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-121.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (127.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (193.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (189.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (79.3892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (63.2539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (152.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (86.3665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (132.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (188.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (131.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-48.9171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-55.1163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-101.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-148.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-227.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-196.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-70.1179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-90.0446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['HF(38)', 'CH2CO(28)', 'C=C[CH]F(1341)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.52349e+10,'s^-1'), n=1.11612, Ea=(203.457,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_N-Sp-7BrCClFINPSSi-2C',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_N-Sp-7BrCClFINPSSi-2C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2F-CH2(62)', 'C=C(O)[C]F(661)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.67744e-07,'m^3/(mol*s)'), n=3.38172, Ea=(60.5478,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', '[CH2]C(F)C(=C)O(673)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(179.647,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['C=C(O)C(F)C[CH]F(10326)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)C=C(O)CF(10327)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(303.681,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['[CH2]C(F)CC(O)=CF(10328)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(161.074,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(18)', 'C=C(O)C(F)[CH]F(678)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', '[CH]C(F)C(F)C(=C)O(10329)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['O[C]1CCC(F)C1F(10330)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.23e+07,'s^-1'), n=1, Ea=(56.484,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_CsCs_RH_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['[CH2]C1(O)CC(F)C1F(10331)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['[CH2]C(F)C(F)C(C)=O(10332)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.913,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C=CC(F)C(=C)O(10333)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(61.9707,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CHF(56)', 'C=C(O)[CH]F(676)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00252,'m^3/(mol*s)'), n=2.41, Ea=(30.8088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_2R!H->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'C=C(O)C(F)C(=C)F(9278)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(45.0881,'m^3/(mol*s)'), n=1.74121, Ea=(17.5012,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_Sp-6R!H-4CClOS_Ext-1COS-R_7R!H->F',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_Sp-6R!H-4CClOS_Ext-1COS-R_7R!H->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[CH2]C(O)=CC([CH2])F(10334)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', '[CH2][CH]C(F)C(=C)O(10335)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OH(7)', '[CH2]C(F)C(F)[C]=C(7449)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[CH2]C(F)C(F)C(=C)[O](10336)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]F(233)', 'C=C(O)[CH]F(676)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(3.70463,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]O(137)', '[CH2]C(F)[CH]F(252)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -23.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[CH2]C(F)[C](F)C(=C)O(10337)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(6.83523,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', '[CH2][C](F)C(F)C(=C)O(10338)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0.261562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', '[CH]=C(O)C(F)C([CH2])F(10339)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F2(78)', '[CH2]C=CC(=C)O(10340)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(22.7805,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', '[CH2]C=C(F)C(=C)O(10341)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(309.319,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', '[CH2]C(F)=CC(=C)O(10342)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(316.527,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(O)C(F)[C](C)F(10343)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.465e+11,'s^-1'), n=0, Ea=(184.64,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 440 used for R2H_S;C_rad_out_noH;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['C=C(O)[C](F)C(C)F(10344)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.59157e+08,'s^-1'), n=1.13, Ea=(121.336,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['C=C([O])C(F)C(C)F(10345)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(O)C(F)C(C)F(10346)'],
    products = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['[CH2]C(O)=CC(F)CF(10347)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(200.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(F)C(F)C(=C)O(10093)'],
    products = ['C=C(O)C(F)[CH]CF(10348)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(180.176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3522',
    isomers = [
        '[CH2]C(F)C(F)C(=C)O(10093)',
    ],
    reactants = [
        ('HF(38)', 'CH2CO(28)', 'C=C[CH]F(1341)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3522',
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

