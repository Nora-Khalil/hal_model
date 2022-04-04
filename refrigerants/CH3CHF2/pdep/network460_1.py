species(
    label = 'O=C(O)C(F)CF(1422)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {11,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-828.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,206,363,518,1175,1317,1481,3059,528,1116,1182,1331,1402,1494,3075,3110,379.296,379.343,379.406,379.421,379.445,2313.33],'cm^-1')),
        HinderedRotor(inertia=(0.0396896,'amu*angstrom^2'), symmetry=1, barrier=(34.6577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339334,'amu*angstrom^2'), symmetry=1, barrier=(34.6577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339462,'amu*angstrom^2'), symmetry=1, barrier=(34.6588,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83821,0.0491707,-3.92463e-05,1.57975e-08,-2.64537e-12,-99617.3,20.2454], Tmin=(100,'K'), Tmax=(1360.04,'K')), NASAPolynomial(coeffs=[10.3348,0.024181,-1.16848e-05,2.28718e-09,-1.61915e-13,-101928,-23.3587], Tmin=(1360.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-828.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs)"""),
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
    label = 'OC(F)CF(544)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {9,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-662.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,528,1116,1182,1331,1402,1494,3075,3110,180],'cm^-1')),
        HinderedRotor(inertia=(1.27049,'amu*angstrom^2'), symmetry=1, barrier=(29.211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27447,'amu*angstrom^2'), symmetry=1, barrier=(29.3025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0494,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89252,0.00611482,0.000108385,-2.11378e-07,1.20595e-10,-79630.2,9.64641], Tmin=(10,'K'), Tmax=(599.394,'K')), NASAPolynomial(coeffs=[3.69664,0.03329,-2.43572e-05,8.26298e-09,-1.04496e-12,-80071.4,6.61484], Tmin=(599.394,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-662.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH3F(39)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-246.136,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.0219,'amu')),
        NonlinearRotor(inertia=([3.21737,19.5168,19.5168],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([1097.73,1206.94,1206.96,1493.65,1493.66,1509.77,3011.42,3086.11,3086.13],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (34.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05191,-0.00340134,1.7341e-05,2.49456e-09,-1.5179e-11,-29605.1,3.91759], Tmin=(10,'K'), Tmax=(534.708,'K')), NASAPolynomial(coeffs=[0.487074,0.0143029,-7.17972e-06,1.7168e-09,-1.57925e-13,-29095.7,20.0824], Tmin=(534.708,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-246.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(O)[C]F(1512)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {5,S}
5 C u0 p1 c0 {1,S} {4,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-238.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,262,1290,465.919,467.467,468.251,469.122,469.638],'cm^-1')),
        HinderedRotor(inertia=(0.000763202,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000753446,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93203,0.0255219,-2.4478e-05,1.25725e-08,-2.73828e-12,-28605.1,14.1558], Tmin=(100,'K'), Tmax=(1066.44,'K')), NASAPolynomial(coeffs=[6.35426,0.012686,-6.42397e-06,1.28649e-09,-9.26041e-14,-29335,-2.57468], Tmin=(1066.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdCsOs) + group(CJ2_singlet-FCO)"""),
)

species(
    label = 'H2(8)',
    structure = adjacencyList("""1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-8.60349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3765.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (2.01594,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(496.376,'J/mol'), sigma=(2.8327,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212708,-2.78619e-07,3.40262e-10,-7.7602e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.09,'K')), NASAPolynomial(coeffs=[2.78813,0.000587684,1.5899e-07,-5.52698e-11,4.34281e-15,-596.125,0.112931], Tmin=(1959.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(O)C(F)[C]F(1513)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u0 p1 c0 {2,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-477.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,617,898,1187,310.225,310.259,310.264,310.271,310.327],'cm^-1')),
        HinderedRotor(inertia=(0.0118243,'amu*angstrom^2'), symmetry=1, barrier=(41.5601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.495092,'amu*angstrom^2'), symmetry=1, barrier=(33.7974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607771,'amu*angstrom^2'), symmetry=1, barrier=(41.5583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98222,0.0489136,-5.66485e-05,3.67444e-08,-1.00684e-11,-57311,20.6882], Tmin=(100,'K'), Tmax=(867.305,'K')), NASAPolynomial(coeffs=[7.58,0.0230965,-1.19973e-05,2.42229e-09,-1.74963e-13,-58282,-5.52084], Tmin=(867.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-477.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C(O)CF(1514)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-603.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,244,1102,1254,1379,1468,1475,3037,3083,180,709.49,709.502,709.546,709.55],'cm^-1')),
        HinderedRotor(inertia=(0.132895,'amu*angstrom^2'), symmetry=1, barrier=(47.4688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132869,'amu*angstrom^2'), symmetry=1, barrier=(47.4695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90573,0.0118318,3.12088e-05,-4.83932e-08,1.93566e-11,-72583,9.24689], Tmin=(10,'K'), Tmax=(885.8,'K')), NASAPolynomial(coeffs=[3.76887,0.0260901,-1.60342e-05,4.54659e-09,-4.90945e-13,-73093.9,6.86993], Tmin=(885.8,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-603.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(O)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144066,5.45064e-06,-3.57997e-09,7.56176e-13,50400.6,-0.411759], Tmin=(100,'K'), Tmax=(1442.38,'K')), NASAPolynomial(coeffs=[2.62651,0.00394757,-1.49921e-06,2.54533e-10,-1.62951e-14,50691.7,6.78356], Tmin=(1442.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(O)C(F)F(1515)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-811.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,195,270,1147,1130,1359,1388,1409,3075,302.397,303.118,303.679,304.78,1108.99],'cm^-1')),
        HinderedRotor(inertia=(0.0134643,'amu*angstrom^2'), symmetry=1, barrier=(29.4664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0134358,'amu*angstrom^2'), symmetry=1, barrier=(29.4625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.00261,0.0220897,-2.96049e-06,-9.04426e-09,3.91266e-12,-97592.5,10.9495], Tmin=(10,'K'), Tmax=(1184.8,'K')), NASAPolynomial(coeffs=[10.5537,0.0121968,-5.91178e-06,1.32461e-09,-1.12746e-13,-100003,-25.3872], Tmin=(1184.8,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-811.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(O)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FCCF(225)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-461.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([442,614,1089,1143,1126,1238,1254,1408,1339,1465,1463,1525,3023,3127,3088,3132,180],'cm^-1')),
        HinderedRotor(inertia=(0.533588,'amu*angstrom^2'), symmetry=1, barrier=(17.8105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.05,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2595.78,'J/mol'), sigma=(4.583,'angstroms'), dipoleMoment=(2,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92422,0.00937296,2.43709e-05,-3.4286e-08,1.28696e-11,-55562,8.62109], Tmin=(10,'K'), Tmax=(910.507,'K')), NASAPolynomial(coeffs=[2.74666,0.0232304,-1.27649e-05,3.37977e-09,-3.48575e-13,-55707.5,12.215], Tmin=(910.507,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-461.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""FCCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C=CCF(1442)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {8,S}
5 C u0 p0 c0 {2,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-266.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.459648,'amu*angstrom^2'), symmetry=1, barrier=(14.4466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90028,0.00862533,9.2651e-05,-2.71363e-07,2.47548e-10,-32065.5,10.3159], Tmin=(10,'K'), Tmax=(342.261,'K')), NASAPolynomial(coeffs=[3.02631,0.0262359,-1.69456e-05,5.25431e-09,-6.24003e-13,-32049,12.9624], Tmin=(342.261,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-266.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""ODCDCCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787932,2.90876e-06,-1.47517e-09,2.12837e-13,-30281.6,-0.311363], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81029e-07,1.07243e-10,-5.79388e-15,-29958.6,5.91041], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C=C(F)CF(1516)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {2,S} {4,S} {6,D}
6 C u0 p0 c0 {3,D} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-395.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,145,326,398,834,1303,2120,512.5,787.5,249.818],'cm^-1')),
        HinderedRotor(inertia=(0.00270418,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73053,0.0251493,-3.61899e-06,-1.84764e-08,1.12937e-11,-47556.6,11.4901], Tmin=(10,'K'), Tmax=(772.96,'K')), NASAPolynomial(coeffs=[5.94217,0.0219841,-1.35443e-05,3.94213e-09,-4.39302e-13,-48145.8,-0.210171], Tmin=(772.96,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-395.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""ODCDC(F)CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OCF(191)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {6,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-436.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,548,1085,1183,1302,1466,1520,3060,3119],'cm^-1')),
        HinderedRotor(inertia=(0.193143,'amu*angstrom^2'), symmetry=1, barrier=(16.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0324,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03993,-0.00338798,6.00172e-05,-9.49101e-08,4.84175e-11,-52506.5,7.78779], Tmin=(10,'K'), Tmax=(607.141,'K')), NASAPolynomial(coeffs=[1.98113,0.0176998,-1.06708e-05,3.11934e-09,-3.52169e-13,-52395.2,15.551], Tmin=(607.141,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-436.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=CF(1070)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-172.285,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(60.0011,'amu')),
        NonlinearRotor(inertia=([9.01649,110.348,119.365],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([238.894,460.044,539.472,686.431,1048.64,1208.61,1432.92,2235.18,3235.83],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95581,0.00278106,4.95039e-05,-1.085e-07,7.1303e-11,-20718.2,7.95963], Tmin=(10,'K'), Tmax=(510.864,'K')), NASAPolynomial(coeffs=[3.75632,0.0135037,-8.87766e-06,2.78796e-09,-3.34853e-13,-20817.4,7.61808], Tmin=(510.864,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-172.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ODCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(=CCF)OF(1445)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {7,S} {11,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-377.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,212,931,1104,1251,1325,1474,3102,3155,3010,987.5,1337.5,450,1655,350,440,435,1725,279.499,279.515],'cm^-1')),
        HinderedRotor(inertia=(0.236404,'amu*angstrom^2'), symmetry=1, barrier=(13.0954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236204,'amu*angstrom^2'), symmetry=1, barrier=(13.0957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251591,'amu*angstrom^2'), symmetry=1, barrier=(13.9491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975728,0.0709771,-0.000101071,7.62928e-08,-2.29954e-11,-45297.9,21.7908], Tmin=(100,'K'), Tmax=(812.177,'K')), NASAPolynomial(coeffs=[10.9183,0.0220092,-1.06311e-05,2.05571e-09,-1.43878e-13,-46912.9,-24.1077], Tmin=(812.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs)"""),
)

species(
    label = 'OC(=CF)OCF(1517)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {6,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-684.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,548,1085,1183,1302,1466,1520,3060,3119,350,440,435,1725,194,682,905,1196,1383,3221,231.605,231.713,231.764],'cm^-1')),
        HinderedRotor(inertia=(0.480885,'amu*angstrom^2'), symmetry=1, barrier=(18.3346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48115,'amu*angstrom^2'), symmetry=1, barrier=(18.3378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.480976,'amu*angstrom^2'), symmetry=1, barrier=(18.3387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921203,0.0657113,-7.48696e-05,4.20867e-08,-9.24111e-12,-82205.9,20.7592], Tmin=(100,'K'), Tmax=(1115.41,'K')), NASAPolynomial(coeffs=[15.0881,0.0149069,-6.54778e-06,1.2515e-09,-8.85658e-14,-85366.2,-49.135], Tmin=(1115.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsFHHO) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C([O])C(F)CF(1518)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {7,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-415.619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,528,1116,1182,1331,1402,1494,3075,3110,1380,1390,370,380,2900,435,180,180,1104.33,1109.8],'cm^-1')),
        HinderedRotor(inertia=(0.229631,'amu*angstrom^2'), symmetry=1, barrier=(5.27966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230582,'amu*angstrom^2'), symmetry=1, barrier=(5.30154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47955,0.0661987,-0.000114958,1.17159e-07,-4.59048e-11,-49907.2,24.4978], Tmin=(100,'K'), Tmax=(817.841,'K')), NASAPolynomial(coeffs=[2.1183,0.0389003,-2.05523e-05,4.06124e-09,-2.84741e-13,-49203.2,26.4873], Tmin=(817.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-415.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(O)[C](F)CF(1519)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {11,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {5,S} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-450.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,180,180,1178.82],'cm^-1')),
        HinderedRotor(inertia=(0.208326,'amu*angstrom^2'), symmetry=1, barrier=(4.78982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209632,'amu*angstrom^2'), symmetry=1, barrier=(4.81985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207437,'amu*angstrom^2'), symmetry=1, barrier=(4.76939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19194,0.0726099,-0.000129465,1.26458e-07,-4.71854e-11,-54134.1,27.3745], Tmin=(100,'K'), Tmax=(835.275,'K')), NASAPolynomial(coeffs=[4.30323,0.0340787,-1.78321e-05,3.48716e-09,-2.42134e-13,-53829.5,17.8591], Tmin=(835.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CCOJ) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(O)C(F)[CH]F(1520)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {11,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
7  C u1 p0 c0 {2,S} {5,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-449.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,1380,1390,370,380,2900,435,334,575,1197,1424,3202,180,180,671.522],'cm^-1')),
        HinderedRotor(inertia=(0.174132,'amu*angstrom^2'), symmetry=1, barrier=(4.00363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182773,'amu*angstrom^2'), symmetry=1, barrier=(4.20232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179781,'amu*angstrom^2'), symmetry=1, barrier=(4.13352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969319,0.0761975,-0.000132217,1.23965e-07,-4.49431e-11,-53910.4,26.6901], Tmin=(100,'K'), Tmax=(828.794,'K')), NASAPolynomial(coeffs=[6.43326,0.030959,-1.61935e-05,3.17006e-09,-2.20445e-13,-54168.1,5.26516], Tmin=(828.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-449.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CCOJ) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'O[C](O)C(F)[CH]F(1521)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {11,S}
4  O u0 p2 c0 {6,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u1 p0 c0 {2,S} {5,S} {9,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-469.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,259,529,569,1128,1321,1390,3140,360,370,350,334,575,1197,1424,3202,229.082,229.232],'cm^-1')),
        HinderedRotor(inertia=(0.00320868,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287712,'amu*angstrom^2'), symmetry=1, barrier=(10.7173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00321193,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.342679,'amu*angstrom^2'), symmetry=1, barrier=(12.7651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674972,0.0792749,-0.00012874,1.08613e-07,-3.59819e-11,-56357.1,26.8834], Tmin=(100,'K'), Tmax=(813.371,'K')), NASAPolynomial(coeffs=[10.8994,0.0222829,-1.12618e-05,2.18136e-09,-1.51152e-13,-57798.4,-18.9671], Tmin=(813.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsOsH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Cs_P) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'OC(O)=C(F)CF(1522)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {11,S}
4  O u0 p2 c0 {7,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {3,S} {4,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-706.687,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,212,931,1104,1251,1325,1474,3102,3155,323,467,575,827,1418,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.726344,'amu*angstrom^2'), symmetry=1, barrier=(16.7001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726799,'amu*angstrom^2'), symmetry=1, barrier=(16.7105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726408,'amu*angstrom^2'), symmetry=1, barrier=(16.7016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777108,0.0722104,-9.51635e-05,6.26664e-08,-1.61084e-11,-84879.9,21.0831], Tmin=(100,'K'), Tmax=(957.349,'K')), NASAPolynomial(coeffs=[14.2297,0.0160026,-7.09527e-06,1.33823e-09,-9.31946e-14,-87455.6,-43.2312], Tmin=(957.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-706.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsCsCs)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(O)[CH]CF(1523)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {6,S} {10,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u1 p0 c0 {4,S} {6,S} {9,S}
6  C u0 p0 c0 {2,S} {3,D} {5,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-457.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,3025,407.5,1350,352.5,327.273,327.275,327.276,327.279,919.631,1693.1],'cm^-1')),
        HinderedRotor(inertia=(0.53107,'amu*angstrom^2'), symmetry=1, barrier=(40.3654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00157394,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0198433,'amu*angstrom^2'), symmetry=1, barrier=(40.3654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40618,0.0377434,-2.63005e-05,9.27077e-09,-1.3924e-12,-55023.2,19.8475], Tmin=(100,'K'), Tmax=(1445.73,'K')), NASAPolynomial(coeffs=[7.86404,0.0226427,-1.0633e-05,2.04605e-09,-1.43073e-13,-56601.3,-8.4953], Tmin=(1445.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-457.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(Cds-OdCsOs) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(F)C(=O)O(1524)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {3,D} {4,S}
6  C u1 p0 c0 {4,S} {8,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-445.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,180,233,1109,1254,1325,1339,3290,3000,3100,440,815,1455,1000,297.039,297.043,297.058,297.059,829.367],'cm^-1')),
        HinderedRotor(inertia=(0.00190991,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684381,'amu*angstrom^2'), symmetry=1, barrier=(42.8617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0297572,'amu*angstrom^2'), symmetry=1, barrier=(42.8615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99662,0.0484039,-5.00907e-05,2.94187e-08,-7.4261e-12,-53555.4,18.8314], Tmin=(100,'K'), Tmax=(927.96,'K')), NASAPolynomial(coeffs=[7.36858,0.0252482,-1.26612e-05,2.52881e-09,-1.8184e-13,-54552.4,-6.68364], Tmin=(927.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(Csj(Cs-F1sCOH)(H)(H))"""),
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
    label = 'O=[C]C(F)CF(1525)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-409.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([255,365,576,1163,1267,1448,3060,528,1116,1182,1331,1402,1494,3075,3110,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.449417,'amu*angstrom^2'), symmetry=1, barrier=(10.333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15412,'amu*angstrom^2'), symmetry=1, barrier=(26.5355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59537,0.0422143,-9.93972e-05,2.1208e-07,-1.8699e-10,-49256.1,11.9557], Tmin=(10,'K'), Tmax=(352.005,'K')), NASAPolynomial(coeffs=[3.98096,0.0308069,-2.08481e-05,6.61412e-09,-7.94088e-13,-49239.8,11.1163], Tmin=(352.005,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-409.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OD[C]C(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(=O)C(F)CF(1526)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-603.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,363,518,1175,1317,1481,3059,528,1116,1182,1331,1402,1494,3075,3110,180,180,652.162,652.17,652.172,652.173,652.173],'cm^-1')),
        HinderedRotor(inertia=(0.00551372,'amu*angstrom^2'), symmetry=1, barrier=(1.66416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176011,'amu*angstrom^2'), symmetry=1, barrier=(53.1242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00476,0.0487614,-6.04255e-05,5.47906e-08,-2.24848e-11,-72480.1,19.818], Tmin=(100,'K'), Tmax=(683.266,'K')), NASAPolynomial(coeffs=[3.80795,0.033754,-1.77077e-05,3.5763e-09,-2.57599e-13,-72622.6,12.5657], Tmin=(683.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-603.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CCOJ)"""),
)

species(
    label = 'CH2F(46)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-42.5685,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(33.0141,'amu')),
        NonlinearRotor(inertia=([1.91548,16.2277,17.9803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([576.418,1180.5,1217.62,1485.55,3118.23,3268.88],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.025,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03338,-0.00262849,2.74227e-05,-3.89096e-08,1.85259e-11,-5119.82,5.20374], Tmin=(10,'K'), Tmax=(594.366,'K')), NASAPolynomial(coeffs=[2.59024,0.00857266,-4.60348e-06,1.22743e-09,-1.29255e-13,-4974.57,11.194], Tmin=(594.366,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-42.5685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[CH2]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(O)[CH]F(1527)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {5,S}
5 C u1 p0 c0 {1,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-437.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,235,1215,1347,1486,3221,452.939,452.944,452.944,452.949,452.954],'cm^-1')),
        HinderedRotor(inertia=(0.325083,'amu*angstrom^2'), symmetry=1, barrier=(47.3279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000821634,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91283,0.00519456,8.72792e-05,-1.82058e-07,1.113e-10,-52625.5,10.0943], Tmin=(10,'K'), Tmax=(558.971,'K')), NASAPolynomial(coeffs=[3.89109,0.024878,-1.79443e-05,5.9387e-09,-7.34256e-13,-52928.2,7.45768], Tmin=(558.971,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-437.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(O)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]O(165)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (-192.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1244.76,1684.94],'cm^-1')),
        HinderedRotor(inertia=(0.00103939,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-192.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2F-CHF(66)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-267.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,334,575,1197,1424,3202,700.651],'cm^-1')),
        HinderedRotor(inertia=(0.57132,'amu*angstrom^2'), symmetry=1, barrier=(13.1358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93382,0.0142501,1.77758e-06,-8.62923e-09,3.21423e-12,-32220.6,9.76326], Tmin=(10,'K'), Tmax=(1157.03,'K')), NASAPolynomial(coeffs=[6.23595,0.0134751,-6.53103e-06,1.52429e-09,-1.39119e-13,-33234.2,-3.75694], Tmin=(1157.03,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-267.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""F[CH]CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(O)[C](F)CF(1528)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u1 p0 c0 {2,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {4,D} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-672.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,234,347,1316,1464,435.609,435.793,435.921,435.989,436.122,1772.18],'cm^-1')),
        HinderedRotor(inertia=(0.373529,'amu*angstrom^2'), symmetry=1, barrier=(50.3431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232331,'amu*angstrom^2'), symmetry=1, barrier=(31.3382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319407,'amu*angstrom^2'), symmetry=1, barrier=(43.0744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83865,0.0471033,-3.97137e-05,1.68486e-08,-2.92321e-12,-80808.6,20.7456], Tmin=(100,'K'), Tmax=(1343.23,'K')), NASAPolynomial(coeffs=[10.9317,0.0200252,-9.47522e-06,1.84075e-09,-1.29966e-13,-83251.4,-25.806], Tmin=(1343.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-672.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CsCOCsF1s)"""),
)

species(
    label = 'O=C(O)C(F)[CH]F(1529)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {4,D} {5,S}
7  C u1 p0 c0 {2,S} {5,S} {9,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-632.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,180,233,1109,1254,1325,1339,3290,334,575,1197,1424,3202,421.932,421.972,422.164,422.175,422.264,422.437],'cm^-1')),
        HinderedRotor(inertia=(0.000946275,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000945546,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416364,'amu*angstrom^2'), symmetry=1, barrier=(52.5868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54537,0.0589702,-8.01387e-05,6.28893e-08,-2.05449e-11,-75945.3,22.4716], Tmin=(100,'K'), Tmax=(740.064,'K')), NASAPolynomial(coeffs=[7.649,0.0259764,-1.32572e-05,2.63361e-09,-1.87545e-13,-76848.6,-5.13685], Tmin=(740.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-632.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(Csj(Cs-F1sCOH)(F1s)(H))"""),
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
    label = 'C=CC(=O)O(1530)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 C u0 p0 c0 {3,D} {7,S} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-287.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.216788,'amu*angstrom^2'), symmetry=1, barrier=(4.98437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216296,'amu*angstrom^2'), symmetry=1, barrier=(4.97306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38944,0.0170998,-8.39552e-06,1.63165e-09,-9.99466e-14,-34570.6,12.3568], Tmin=(100,'K'), Tmax=(2312.87,'K')), NASAPolynomial(coeffs=[11.8116,0.0053,-2.53658e-06,4.599e-10,-2.91781e-14,-39206.2,-36.9362], Tmin=(2312.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(F)C(=O)O(176)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,D}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 C u0 p0 c0 {4,D} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-451.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,288,410,724,839,1320,2950,3100,1380,975,1025,1650,312.481,312.492,2228.9,2228.91,2228.93],'cm^-1')),
        HinderedRotor(inertia=(0.43152,'amu*angstrom^2'), symmetry=1, barrier=(29.9045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159697,'amu*angstrom^2'), symmetry=1, barrier=(11.0668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53159,0.0354104,-3.27009e-05,1.82662e-08,-4.57503e-12,-54218.2,17.8793], Tmin=(100,'K'), Tmax=(915.684,'K')), NASAPolynomial(coeffs=[5.64595,0.021806,-1.04155e-05,2.04135e-09,-1.45358e-13,-54788.6,3.12859], Tmin=(915.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-451.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'O=C(O)C=CF(175)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 C u0 p0 c0 {1,S} {4,D} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-479.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.159391,'amu*angstrom^2'), symmetry=1, barrier=(3.66472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157813,'amu*angstrom^2'), symmetry=1, barrier=(3.62843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72066,0.0314381,-4.6626e-05,4.33273e-08,-1.66054e-11,-57613.3,15.6138], Tmin=(100,'K'), Tmax=(747.204,'K')), NASAPolynomial(coeffs=[4.66449,0.0167576,-8.57412e-06,1.72068e-09,-1.23035e-13,-57784.4,7.60092], Tmin=(747.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFH)"""),
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
    E0 = (-324.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-260.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-129.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-100.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-138.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (63.9999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-239.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (113.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-174.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-56.0082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (55.6518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-283.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-89.0685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-124.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-81.9779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-131.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-238.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-78.3721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-69.2733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-77.4549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-87.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-176.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-156.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-149.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-116.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (13.5752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-228.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(O)C(F)CF(1422)'],
    products = ['HF(38)', 'CO2(14)', 'CH2CHF(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.30828e+07,'s^-1'), n=1.83354, Ea=(200.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.49335582114639515, var=20.078174816455903, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_7R!H->O',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_7R!H->O"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'OC(F)CF(544)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(216.721,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH3F(39)', 'O=C(O)[C]F(1512)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.00162e-06,'m^3/(mol*s)'), n=3.38172, Ea=(50.7513,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H2(8)', 'O=C(O)C(F)[C]F(1513)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.26413e-09,'m^3/(mol*s)'), n=4.30786, Ea=(81.6268,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HH_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node HH_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHF(40)', 'O=C(O)CF(1514)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.59914e-05,'m^3/(mol*s)'), n=3.10993, Ea=(22.5157,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_N-2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(S)(25)', 'O=C(O)C(F)F(1515)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(152.586,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CO2(14)', 'FCCF(225)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(848000,'cm^3/(mol*s)'), n=2.13, Ea=(322.168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [CO2_Cdd;C_sec]
Euclidian distance = 0
Multiplied by reaction path degeneracy 8.0
family: 1,3_Insertion_CO2"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OF(142)', 'O=C=CCF(1442)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(51.5,'cm^3/(mol*s)'), n=3.05, Ea=(171.544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for cco_HNd;H_OH
Exact match found for rate rule [cco_HNd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2O(3)', 'O=C=C(F)CF(1516)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.000454348,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_Cdd;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OCF(191)', 'O=C=CF(1070)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.37684e-08,'m^3/(mol*s)'), n=3.46667, Ea=(249.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [doublebond;R_OH] + [Cd_Cdd;R_OR] for rate rule [Cd_Cdd;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC(=CCF)OF(1445)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(129.462,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OC(=CF)OCF(1517)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(97.3511,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C([O])C(F)CF(1518)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(O)[C](F)CF(1519)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(O)C(F)[CH]F(1520)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O[C](O)C(F)[CH]F(1521)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OC(O)=C(F)CF(1522)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(410000,'s^-1'), n=2.37, Ea=(164.067,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', 'O=C(O)[CH]CF(1523)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(2.99009,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F(37)', '[CH2]C(F)C(=O)O(1524)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OH(7)', 'O=[C]C(F)CF(1525)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]C(=O)C(F)CF(1526)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2F(46)', 'O=C(O)[CH]F(1527)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -2.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]O(165)', 'CH2F-CHF(66)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', 'O=C(O)[C](F)CF(1528)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.57075e+24,'m^3/(mol*s)'), n=-6.76739, Ea=(7.20789,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-2C-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(5)', 'O=C(O)C(F)[CH]F(1529)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['F2(78)', 'C=CC(=O)O(1530)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(6.28045,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', 'C=C(F)C(=O)O(176)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(52.9173,'m^3/(mol*s)'), n=1.09798, Ea=(200.597,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_3COCdCddCtO2d->Cd_Ext-3Cd-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_3COCdCddCtO2d->Cd_Ext-3Cd-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['HF(38)', 'O=C(O)C=CF(175)'],
    products = ['O=C(O)C(F)CF(1422)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(230.804,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #460',
    isomers = [
        'O=C(O)C(F)CF(1422)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'CH2CHF(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #460',
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

