species(
    label = '[O]CC(F)C(=O)OF(10805)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-454.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,206,363,518,1175,1317,1481,3059,2750,2850,1437.5,1250,1305,750,350,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13238,0.0718154,-0.000116184,1.11191e-07,-4.24821e-11,-54550.9,25.2039], Tmin=(100,'K'), Tmax=(786.647,'K')), NASAPolynomial(coeffs=[4.85546,0.037073,-1.97874e-05,3.94683e-09,-2.79263e-13,-54647.5,11.2451], Tmin=(786.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-454.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ)"""),
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
    label = 'CO2(14)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C[CH]F(338)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-191.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.21522,'amu*angstrom^2'), symmetry=1, barrier=(35.2035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2974.27,'J/mol'), sigma=(4.91649,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.57 K, Pc=56.79 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96722,0.00209998,6.04576e-05,-1.2409e-07,8.01311e-11,-23018.4,8.77531], Tmin=(10,'K'), Tmax=(478.869,'K')), NASAPolynomial(coeffs=[2.63637,0.0192853,-1.23829e-05,3.78104e-09,-4.41625e-13,-22960.5,13.4894], Tmin=(478.869,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-191.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CO(13)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
"""),
    E0 = (-118.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2193.04],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]CC(F)OF(1745)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-283.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,261,493,600,1152,1365,1422,3097,2750,2850,1437.5,1250,1305,750,350,212.833,213.001,1715.25],'cm^-1')),
        HinderedRotor(inertia=(0.721553,'amu*angstrom^2'), symmetry=1, barrier=(22.9465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01011,'amu*angstrom^2'), symmetry=1, barrier=(32.2246,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60793,0.0359301,-2.69116e-05,8.87659e-09,-8.98365e-13,-34087,12.2367], Tmin=(10,'K'), Tmax=(1178.51,'K')), NASAPolynomial(coeffs=[10.0899,0.0184081,-1.031e-05,2.70989e-09,-2.74257e-13,-35925.8,-21.4192], Tmin=(1178.51,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-283.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]CC(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH3O(27)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (10.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.15,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71181,-0.00280463,3.76551e-05,-4.73072e-08,1.86588e-11,1295.7,6.57241], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.75779,0.00744142,-2.69705e-06,4.38091e-10,-2.63537e-14,378.112,-1.9668], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(10.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH3O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C([C]F)OF(1423)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {3,S} {4,D} {6,S}
6 C u0 p1 c0 {1,S} {5,S}
"""),
    E0 = (-103.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,262,1290,180,180,675.255,1002.35,1005.96,1006.9],'cm^-1')),
        HinderedRotor(inertia=(2.3037,'amu*angstrom^2'), symmetry=1, barrier=(52.9665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30448,'amu*angstrom^2'), symmetry=1, barrier=(52.9845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44534,0.0378912,-5.06822e-05,3.61209e-08,-1.05858e-11,-12411.6,16.0814], Tmin=(100,'K'), Tmax=(822.606,'K')), NASAPolynomial(coeffs=[7.23663,0.0145939,-8.2018e-06,1.69455e-09,-1.23553e-13,-13199.8,-6.09837], Tmin=(822.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cds-OdCsOs) + group(CJ2_singlet-FCO)"""),
)

species(
    label = 'FOF(499)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {1,S} {2,S}
"""),
    E0 = (16.0257,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(53.9917,'amu')),
        NonlinearRotor(inertia=([8.34135,45.8552,54.1965],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([486.983,915.713,1034.58],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (53.9962,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03449,-0.00336126,4.20188e-05,-7.82672e-08,4.57369e-11,1928.28,6.37844], Tmin=(10,'K'), Tmax=(574.663,'K')), NASAPolynomial(coeffs=[3.91346,0.00598694,-4.58407e-06,1.5533e-09,-1.93094e-13,1801.75,5.67331], Tmin=(574.663,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(16.0257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""FOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]CC=C=O(913)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {8,S}
5 C u0 p0 c0 {2,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-11.8211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0898161,'amu*angstrom^2'), symmetry=1, barrier=(2.06505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43868,0.0419871,-7.44652e-05,7.48279e-08,-2.80499e-11,-1372.86,16.4118], Tmin=(100,'K'), Tmax=(869.697,'K')), NASAPolynomial(coeffs=[2.35089,0.023934,-1.14954e-05,2.15667e-09,-1.4569e-13,-659.577,20.8361], Tmin=(869.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.8211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CCOJ)"""),
)

species(
    label = 'OF(142)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]CC(F)=C=O(5428)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u0 p0 c0 {3,D} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-161.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,145,326,398,834,1303,2120,512.5,787.5,180,2321.18],'cm^-1')),
        HinderedRotor(inertia=(0.346485,'amu*angstrom^2'), symmetry=1, barrier=(7.96638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0451,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92128,0.0554324,-0.000105876,1.04037e-07,-3.78596e-11,-19385.9,18.8212], Tmin=(100,'K'), Tmax=(874.879,'K')), NASAPolynomial(coeffs=[3.81095,0.0237015,-1.18824e-05,2.24458e-09,-1.51244e-13,-18832.9,15.0078], Tmin=(874.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cd(Cdd-Od)CF) + missing(Cdd-CdO2d) + radical(CCOJ)"""),
)

species(
    label = '[O]CC=C(OF)OF(10895)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {1,S} {8,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (28.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([175,287,726,856,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,1678.54,1682.31],'cm^-1')),
        HinderedRotor(inertia=(0.195434,'amu*angstrom^2'), symmetry=1, barrier=(4.49341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196491,'amu*angstrom^2'), symmetry=1, barrier=(4.51772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197101,'amu*angstrom^2'), symmetry=1, barrier=(4.53174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852994,0.0862133,-0.000172599,1.75812e-07,-6.58819e-11,3509.84,27.6528], Tmin=(100,'K'), Tmax=(862.427,'K')), NASAPolynomial(coeffs=[2.14207,0.0396473,-2.10158e-05,4.06783e-09,-2.78336e-13,4796.89,30.3754], Tmin=(862.427,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2sCF) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(CCOJ)"""),
)

species(
    label = '[O]COC(=CF)OF(10896)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {1,S} {7,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-276.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,194,682,905,1196,1383,3221,442.094,442.142,442.168,442.217,1589.34],'cm^-1')),
        HinderedRotor(inertia=(0.0662928,'amu*angstrom^2'), symmetry=1, barrier=(9.19483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200066,'amu*angstrom^2'), symmetry=1, barrier=(27.753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200071,'amu*angstrom^2'), symmetry=1, barrier=(27.7529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0076,0.0707666,-8.67593e-05,5.50746e-08,-1.41561e-11,-33193.2,23.7126], Tmin=(100,'K'), Tmax=(939.05,'K')), NASAPolynomial(coeffs=[11.869,0.024502,-1.28597e-05,2.61164e-09,-1.89365e-13,-35233.2,-28.0046], Tmin=(939.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(OCOJ)"""),
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
    label = '[CH2]C(F)C(=O)OF(1438)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {4,D} {5,S}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-311.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,180,233,1109,1254,1325,1339,3290,3000,3100,440,815,1455,1000,376.569,377.122,377.673,377.744,378.348,378.435],'cm^-1')),
        HinderedRotor(inertia=(0.472592,'amu*angstrom^2'), symmetry=1, barrier=(47.7033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00117412,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47155,'amu*angstrom^2'), symmetry=1, barrier=(47.6509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48006,0.0611777,-7.79421e-05,5.54195e-08,-1.64571e-11,-37360.7,20.8604], Tmin=(100,'K'), Tmax=(808.094,'K')), NASAPolynomial(coeffs=[8.41078,0.0268739,-1.42719e-05,2.89664e-09,-2.09404e-13,-38480.9,-11.1001], Tmin=(808.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(Csj(Cs-F1sCOH)(H)(H))"""),
)

species(
    label = 'FO[C]1OOCC1F(10897)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {2,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-197.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374037,0.0640417,-6.40715e-05,3.18204e-08,-5.93057e-12,-23571.4,25.193], Tmin=(100,'K'), Tmax=(1523.29,'K')), NASAPolynomial(coeffs=[16.8707,0.0105736,-1.42626e-06,2.95877e-11,4.78178e-15,-27419.7,-57.4714], Tmin=(1523.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2sCF) + group(CsCsCsFH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxolane) + radical(Cs_P)"""),
)

species(
    label = '[O]C1(OF)OCC1F(10898)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-311.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18522,0.0462073,-3.18665e-05,-4.48487e-09,1.6913e-11,-37399.2,22.5013], Tmin=(100,'K'), Tmax=(510.437,'K')), NASAPolynomial(coeffs=[4.6269,0.034139,-1.71653e-05,3.43296e-09,-2.47082e-13,-37740.5,11.462], Tmin=(510.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(CsCsCsFH) + group(Cs-CsOsHH) + ring(Cs-O2s-Cs-Cs(F)) + radical(CCOJ)"""),
)

species(
    label = '[O]CC(F)=C(O)OF(10899)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,D}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-311.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,350,440,435,1725,180,180,2044.3],'cm^-1')),
        HinderedRotor(inertia=(0.210059,'amu*angstrom^2'), symmetry=1, barrier=(4.82967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210227,'amu*angstrom^2'), symmetry=1, barrier=(4.83353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210217,'amu*angstrom^2'), symmetry=1, barrier=(4.83331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.618692,0.089368,-0.000173383,1.70017e-07,-6.20805e-11,-37313.5,26.732], Tmin=(100,'K'), Tmax=(861.123,'K')), NASAPolynomial(coeffs=[4.88676,0.0350141,-1.85575e-05,3.59084e-09,-2.45632e-13,-36768.4,14.2125], Tmin=(861.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsCsCs) + radical(CCOJ)"""),
)

species(
    label = 'CH2O(20)',
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
    label = 'O=C([CH]F)OF(1436)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {3,S} {4,D} {6,S}
6 C u1 p0 c0 {1,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-273.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,235,1215,1347,1486,3221,253.617,253.644,253.691,253.72,593.372,2118.11],'cm^-1')),
        HinderedRotor(inertia=(0.691458,'amu*angstrom^2'), symmetry=1, barrier=(31.5756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691582,'amu*angstrom^2'), symmetry=1, barrier=(31.5756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.72467,0.023965,8.45138e-05,-3.95579e-07,4.39047e-10,-32869.5,11.3176], Tmin=(10,'K'), Tmax=(343.968,'K')), NASAPolynomial(coeffs=[5.85162,0.0233125,-1.76582e-05,5.98709e-09,-7.51185e-13,-33158.3,1.25533], Tmin=(343.968,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-273.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""ODC([CH]F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=CC(F)C(=O)OF(10900)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,D} {6,S}
8  C u0 p0 c0 {5,D} {6,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-589.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,207.664,207.768,537.288,537.351,537.443,537.446],'cm^-1')),
        HinderedRotor(inertia=(1.36304,'amu*angstrom^2'), symmetry=1, barrier=(41.747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36076,'amu*angstrom^2'), symmetry=1, barrier=(41.7469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36381,'amu*angstrom^2'), symmetry=1, barrier=(41.7471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61252,0.0583866,-6.80969e-05,4.44654e-08,-1.23243e-11,-70850.4,22.1059], Tmin=(100,'K'), Tmax=(855.567,'K')), NASAPolynomial(coeffs=[8.0976,0.0280671,-1.494e-05,3.045e-09,-2.21036e-13,-71960.1,-8.16932], Tmin=(855.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-589.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
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
    label = '[O]C[CH]C(=O)OF(2729)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u1 p0 c0 {5,S} {7,S} {10,S}
7  C u0 p0 c0 {2,S} {4,D} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-68.7221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74023,0.0575732,-9.23213e-05,9.08168e-08,-3.53344e-11,-8191.66,25.038], Tmin=(100,'K'), Tmax=(802.402,'K')), NASAPolynomial(coeffs=[3.22259,0.0339027,-1.76364e-05,3.47882e-09,-2.44458e-13,-7905.43,21.4788], Tmin=(802.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.7221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = '[O]CC(F)C([O])=O(10901)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u1 p2 c0 {6,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-363.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,363,518,1175,1317,1481,3059,2750,2850,1437.5,1250,1305,750,350,180,180,1011.19,1011.38,1011.39,1011.43,1011.44,1011.49],'cm^-1')),
        HinderedRotor(inertia=(0.130141,'amu*angstrom^2'), symmetry=1, barrier=(2.99219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130171,'amu*angstrom^2'), symmetry=1, barrier=(2.99289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94135,0.0570738,-0.000103115,1.13126e-07,-4.65766e-11,-43614.5,22.2951], Tmin=(100,'K'), Tmax=(819.94,'K')), NASAPolynomial(coeffs=[-0.87839,0.0418872,-2.23848e-05,4.43777e-09,-3.1157e-13,-42179.2,41.2717], Tmin=(819.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]F(126)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]CC(F)[C]=O(10902)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-169.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([255,365,576,1163,1267,1448,3060,2750,2850,1437.5,1250,1305,750,350,1855,455,950,180,1622.66],'cm^-1')),
        HinderedRotor(inertia=(0.241286,'amu*angstrom^2'), symmetry=1, barrier=(5.54764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241122,'amu*angstrom^2'), symmetry=1, barrier=(5.54387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05003,0.0528887,-9.76123e-05,1.01924e-07,-3.99635e-11,-20273.6,20.6198], Tmin=(100,'K'), Tmax=(834.386,'K')), NASAPolynomial(coeffs=[1.57785,0.0315041,-1.66554e-05,3.27273e-09,-2.2794e-13,-19371.6,27.7453], Tmin=(834.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(COj(Cs-F1sCsH)(O2d))"""),
)

species(
    label = '[CH2][O](132)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88412,-0.00363928,3.2856e-05,-4.13634e-08,1.59641e-11,23210.8,7.4797], Tmin=(100,'K'), Tmax=(933.047,'K')), NASAPolynomial(coeffs=[6.69324,0.000290178,8.61305e-07,-1.56324e-10,7.33557e-15,21991.3,-9.60369], Tmin=(933.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = 'O=[C]OF(167)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {2,S} {3,D}
"""),
    E0 = (-35.6837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.77736,'amu*angstrom^2'), symmetry=1, barrier=(40.8651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98172,0.0254564,-4.70452e-05,4.25377e-08,-1.44427e-11,-4258,9.5319], Tmin=(100,'K'), Tmax=(878.471,'K')), NASAPolynomial(coeffs=[5.67748,0.00648671,-3.22256e-06,6.05542e-10,-4.04968e-14,-4473.3,-1.65399], Tmin=(878.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.6837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical((O)CJOC)"""),
)

species(
    label = '[O]C[CH]F(140)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-16.6263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,180,898.208],'cm^-1')),
        HinderedRotor(inertia=(0.176227,'amu*angstrom^2'), symmetry=1, barrier=(14.6052,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3194.99,'J/mol'), sigma=(5.39809,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=499.05 K, Pc=46.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.9265,0.0263719,-2.23156e-05,9.9703e-09,-1.9036e-12,-1963.41,13.4259], Tmin=(100,'K'), Tmax=(1183.98,'K')), NASAPolynomial(coeffs=[6.46068,0.0144322,-7.18937e-06,1.45336e-09,-1.05268e-13,-2800.31,-4.22144], Tmin=(1183.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.6263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = '[O]C[C](F)C(=O)OF(5434)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {1,S} {6,S} {8,S}
8  C u0 p0 c0 {3,S} {5,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-297.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2850,1437.5,1250,1305,750,350,234,347,1316,1464,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26254,0.0682281,-0.000111282,1.04991e-07,-3.94375e-11,-35748,25.2371], Tmin=(100,'K'), Tmax=(790.924,'K')), NASAPolynomial(coeffs=[5.47196,0.0328711,-1.75455e-05,3.49171e-09,-2.46526e-13,-35973.9,8.69817], Tmin=(790.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-297.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CsCOCsF1s)"""),
)

species(
    label = '[O][CH]C(F)C(=O)OF(10903)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {7,D}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,D} {6,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-274.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,180,233,1109,1254,1325,1339,3290,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861023,0.0805541,-0.000147181,1.4269e-07,-5.29705e-11,-32859,25.7279], Tmin=(100,'K'), Tmax=(826.041,'K')), NASAPolynomial(coeffs=[5.58926,0.0336704,-1.84862e-05,3.67085e-09,-2.56875e-13,-32821.7,8.7742], Tmin=(826.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-274.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C=CC(=O)OF(2735)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {7,S}
5 C u0 p0 c0 {6,S} {7,D} {8,S}
6 C u0 p0 c0 {2,S} {3,D} {5,S}
7 C u0 p0 c0 {4,S} {5,D} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-220.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.25088,'amu*angstrom^2'), symmetry=1, barrier=(5.76822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249541,'amu*angstrom^2'), symmetry=1, barrier=(5.73743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03357,0.0411722,-4.61732e-05,2.43096e-08,-4.93303e-12,-26436.8,18.6703], Tmin=(100,'K'), Tmax=(1208.87,'K')), NASAPolynomial(coeffs=[12.4535,0.00669359,-3.39087e-06,7.15869e-10,-5.36919e-14,-28956.1,-33.5764], Tmin=(1208.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'O=C(OF)C(F)[CH]O(10904)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u1 p0 c0 {4,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-499.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3615,1277.5,1000,180,233,1109,1254,1325,1339,3290,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.565359,0.0824572,-0.000131658,1.13022e-07,-3.8742e-11,-59989.8,25.9326], Tmin=(100,'K'), Tmax=(777.166,'K')), NASAPolynomial(coeffs=[10.0931,0.0275879,-1.45009e-05,2.86777e-09,-2.01819e-13,-61294.6,-16.4982], Tmin=(777.166,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-499.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCsJOH)"""),
)

species(
    label = 'O=C(OF)[C](F)CO(10905)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {11,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {1,S} {6,S} {8,S}
8  C u0 p0 c0 {4,S} {5,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-523.686,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,2750,2850,1437.5,1250,1305,750,350,234,347,1316,1464,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17729,0.067372,-8.46232e-05,5.83049e-08,-1.65491e-11,-62887.8,24.7039], Tmin=(100,'K'), Tmax=(849.119,'K')), NASAPolynomial(coeffs=[9.67829,0.027327,-1.38845e-05,2.76767e-09,-1.98178e-13,-64331.6,-14.9184], Tmin=(849.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-523.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CsCOCsF1s)"""),
)

species(
    label = '[O]C(=O)C(F)COF(10906)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-454.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,206,363,518,1175,1317,1481,3059,2750,2850,1437.5,1250,1305,750,350,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13238,0.0718154,-0.000116184,1.11191e-07,-4.24821e-11,-54550.9,24.5108], Tmin=(100,'K'), Tmax=(786.647,'K')), NASAPolynomial(coeffs=[4.85546,0.037073,-1.97874e-05,3.94683e-09,-2.79263e-13,-54647.5,10.5519], Tmin=(786.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-454.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ)"""),
)

species(
    label = 'O=C([CH]COF)OF(10907)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {1,S} {6,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {6,S} {8,S} {11,S}
8  C u0 p0 c0 {4,S} {5,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-159.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,990,1113,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07214,0.0704225,-9.75471e-05,7.66354e-08,-2.49143e-11,-19134,26.7624], Tmin=(100,'K'), Tmax=(744.759,'K')), NASAPolynomial(coeffs=[8.66975,0.0296167,-1.53612e-05,3.06715e-09,-2.18927e-13,-20265.7,-7.65276], Tmin=(744.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2sCF) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCJCO)"""),
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
    E0 = (-286.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-30.7142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (111.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (332.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (69.1466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (306.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-33.7138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (88.7952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-40.0741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-154.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1.74844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-216.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-203.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (161.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-133.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (90.7264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (76.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (104.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (70.9313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (94.8603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-56.4315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-213.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-297.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-99.0213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (104.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CC(F)C(=O)OF(10805)'],
    products = ['HF(38)', 'CO2(14)', 'O=C[CH]F(338)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(11.2412,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]CC(F)OF(1745)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(214.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH3O(27)', 'O=C([C]F)OF(1423)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.00162e-06,'m^3/(mol*s)'), n=3.38172, Ea=(47.1132,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['FOF(499)', '[O]CC=C=O(913)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(103,'cm^3/(mol*s)'), n=3.05, Ea=(171.544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for cco_HNd;H_OH
Exact match found for rate rule [cco_HNd;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(142)', '[O]CC(F)=C=O(5428)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000227174,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_Cdd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]CC=C(OF)OF(10895)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(120.796,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]COC(=CF)OF(10896)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(86.0245,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(6)', '[CH2]C(F)C(=O)OF(1438)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]CC(F)C(=O)OF(10805)'],
    products = ['FO[C]1OOCC1F(10897)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(257.168,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 255.5 to 257.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]CC(F)C(=O)OF(10805)'],
    products = ['[O]C1(OF)OCC1F(10898)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(142.903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 140.3 to 142.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]CC(F)=C(O)OF(10899)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(155.776,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2O(20)', 'O=C([CH]F)OF(1436)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.22525e-05,'m^3/(mol*s)'), n=2.9005, Ea=(18.8849,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_N-2R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'O=CC(F)C(=O)OF(10900)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.99,'m^3/(mol*s)'), n=2.12, Ea=(17.3626,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_N-Sp-5R!H=4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_N-Sp-5R!H=4C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]C[CH]C(=O)OF(2729)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]CC(F)C([O])=O(10901)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.62e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]F(126)', '[O]CC(F)[C]=O(10902)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][O](132)', 'O=C([CH]F)OF(1436)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]OF(167)', '[O]C[CH]F(140)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[O]C[C](F)C(=O)OF(5434)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O][CH]C(F)C(=O)OF(10903)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C=CC(=O)OF(2735)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(287.986,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C(OF)C(F)[CH]O(10904)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C(OF)[C](F)CO(10905)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.248e+06,'s^-1'), n=0.992, Ea=(69.3363,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 441 used for R3H_SS_Cs;C_rad_out_noH;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_noH;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 69.3 to 69.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)C(F)COF(10906)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(198.221,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C([CH]COF)OF(10907)'],
    products = ['[O]CC(F)C(=O)OF(10805)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(107.239,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #3775',
    isomers = [
        '[O]CC(F)C(=O)OF(10805)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'O=C[CH]F(338)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3775',
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

