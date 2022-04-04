species(
    label = '[CH]=C(CF)C(F)[CH]F(9331)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {4,S} {12,S}
8  C u1 p0 c0 {6,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-164.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,334,575,1197,1424,3202,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.21437,'amu*angstrom^2'), symmetry=1, barrier=(27.9208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418543,'amu*angstrom^2'), symmetry=1, barrier=(9.62313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.4189,'amu*angstrom^2'), symmetry=1, barrier=(9.63134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40802,0.0857582,-0.000128003,1.06951e-07,-3.62914e-11,-19636.6,27.7018], Tmin=(100,'K'), Tmax=(755.904,'K')), NASAPolynomial(coeffs=[9.88779,0.0325582,-1.64092e-05,3.21809e-09,-2.26361e-13,-20983.1,-14.8061], Tmin=(755.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(Cds_P)"""),
)

species(
    label = 'CHFCHF[Z](58)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-310.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([151,237,609,755,844,966,1147,1245,1323,1443,3181,3261],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383152,0.0294896,-2.94145e-05,1.64336e-08,-4.01759e-12,-36926.9,22.5083], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.34201,0.00821939,-3.17549e-06,5.49282e-10,-3.47434e-14,-38823.3,-13.1129], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-310.115,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCHF[Z]""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'C#CCF(4379)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (8.12032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([319,1023,1071,1259,1317,1409,3054,3019,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.60388,'amu*angstrom^2'), symmetry=1, barrier=(36.8763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2742,'J/mol'), sigma=(4.732,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.29 K, Pc=58.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95512,0.00280826,7.03906e-05,-1.43186e-07,9.06266e-11,979.285,8.14914], Tmin=(10,'K'), Tmax=(503.334,'K')), NASAPolynomial(coeffs=[2.71827,0.0217978,-1.34993e-05,4.0832e-09,-4.79139e-13,987.76,12.1145], Tmin=(503.334,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(8.12032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""C#CCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]=CC(F)[CH]F(5593)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u1 p0 c0 {2,S} {3,S} {9,S}
6  C u1 p0 c0 {4,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (38.7599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.644532,'amu*angstrom^2'), symmetry=1, barrier=(14.8191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645235,'amu*angstrom^2'), symmetry=1, barrier=(14.8352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0713,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2913.02,'J/mol'), sigma=(5.02618,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=455.01 K, Pc=52.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43693,0.0617059,-9.58854e-05,8.2678e-08,-2.82651e-11,4748.91,17.4825], Tmin=(100,'K'), Tmax=(815.966,'K')), NASAPolynomial(coeffs=[7.81284,0.0229878,-1.09912e-05,2.109e-09,-1.45987e-13,3956.82,-10.4586], Tmin=(815.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.7599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(Cds_P)"""),
)

species(
    label = 'CH2(S)(25)',
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
    label = '[CH]=C(F)C(F)[CH]F(5594)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {4,S} {9,S}
7  C u1 p0 c0 {5,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-140.747,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,246,474,533,1155,334,575,1197,1424,3202,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.24737,'amu*angstrom^2'), symmetry=1, barrier=(5.68752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.804513,'amu*angstrom^2'), symmetry=1, barrier=(18.4973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2931.57,'J/mol'), sigma=(4.93944,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.90 K, Pc=55.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931207,0.0745327,-0.000124783,1.09884e-07,-3.77738e-11,-16824.2,23.0629], Tmin=(100,'K'), Tmax=(820.305,'K')), NASAPolynomial(coeffs=[9.20543,0.0235117,-1.19677e-05,2.33637e-09,-1.6255e-13,-17822.5,-13.0275], Tmin=(820.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]F(837)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.000263365,8.89188e-06,-1.03032e-08,3.5081e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164132,-7.73117e-07,1.90988e-10,-1.59926e-14,25413.4,-0.815515], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C(=CF)CF(11163)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {4,D} {9,S}
6  C u2 p0 c0 {4,S} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-34.8556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,194,682,905,1196,1383,3221,516.27,516.911,518.507,520.795],'cm^-1')),
        HinderedRotor(inertia=(0.270993,'amu*angstrom^2'), symmetry=1, barrier=(52.148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270303,'amu*angstrom^2'), symmetry=1, barrier=(52.1854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0713,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48382,0.0395379,-2.24072e-05,5.97211e-09,-6.57982e-13,-4144.22,17.8821], Tmin=(100,'K'), Tmax=(1831.88,'K')), NASAPolynomial(coeffs=[8.13932,0.0271889,-1.22955e-05,2.29221e-09,-1.55783e-13,-6216.26,-12.826], Tmin=(1831.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.8556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = 'FCC1=CC(F)C1F(11304)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {4,S} {6,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-453.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51234,0.052722,-3.45949e-05,1.05683e-08,-1.28734e-12,-54458.6,21.3099], Tmin=(100,'K'), Tmax=(1856.67,'K')), NASAPolynomial(coeffs=[15.1647,0.0233093,-1.08324e-05,2.03596e-09,-1.38454e-13,-59528.2,-53.0024], Tmin=(1856.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-453.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs(F)-Cs-Cd-Cd)"""),
)

species(
    label = 'C=C(CF)C(F)=CF(9332)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {6,D} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-481.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05819,0.0681051,-6.83095e-05,3.63012e-08,-7.97658e-12,-57763.6,22.1544], Tmin=(100,'K'), Tmax=(1076.75,'K')), NASAPolynomial(coeffs=[11.5094,0.0292802,-1.42235e-05,2.8141e-09,-2.01566e-13,-60014.3,-29.0394], Tmin=(1076.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-481.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
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
    label = '[CH]=C(C=CF)CF(11142)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {2,S} {5,D} {11,S}
7  C u1 p0 c0 {4,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-63.8645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.899093,'amu*angstrom^2'), symmetry=1, barrier=(20.6719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.897373,'amu*angstrom^2'), symmetry=1, barrier=(20.6324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0319,0.0617593,-6.06178e-05,3.03236e-08,-6.02638e-12,-7571.3,22.1569], Tmin=(100,'K'), Tmax=(1217.49,'K')), NASAPolynomial(coeffs=[14.0481,0.0189955,-7.93115e-06,1.47388e-09,-1.02379e-13,-10740.7,-43.2002], Tmin=(1217.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.8645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(CdCFH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CF(4380)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u1 p0 c0 {2,S} {4,D}
4 C u1 p0 c0 {3,D} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (311.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,377.76,379.216,382.52,1900.36,1902.09,1902.61,1902.73,3346.48],'cm^-1')),
        HinderedRotor(inertia=(0.11674,'amu*angstrom^2'), symmetry=1, barrier=(12.0012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74483,0.0317272,-5.02124e-05,4.76765e-08,-1.76459e-11,37561.4,13.6835], Tmin=(100,'K'), Tmax=(834.526,'K')), NASAPolynomial(coeffs=[4.0464,0.0168537,-7.95764e-06,1.52224e-09,-1.05087e-13,37644.8,9.44116], Tmin=(834.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=C(CF)C(F)=CF(11305)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {8,D}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {11,S}
8  C u1 p0 c0 {5,D} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-234.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,280,518,736,852,873,194,682,905,1196,1383,3221,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.703791,'amu*angstrom^2'), symmetry=1, barrier=(16.1815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703625,'amu*angstrom^2'), symmetry=1, barrier=(16.1777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940748,0.072041,-8.59874e-05,5.45295e-08,-1.4099e-11,-28041.8,22.9863], Tmin=(100,'K'), Tmax=(931.807,'K')), NASAPolynomial(coeffs=[11.4226,0.0270443,-1.35518e-05,2.70415e-09,-1.94275e-13,-29995.2,-26.8422], Tmin=(931.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-234.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cds_P)"""),
)

species(
    label = 'F[CH][CH]F(4261)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-71.0739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,406,528,622,1148,1246,1368,1480,3164,3240,1663.85],'cm^-1')),
        HinderedRotor(inertia=(0.367109,'amu*angstrom^2'), symmetry=1, barrier=(8.44055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88778,0.0292785,-5.20594e-05,5.24809e-08,-2.00716e-11,-8512.84,14.1703], Tmin=(100,'K'), Tmax=(832.93,'K')), NASAPolynomial(coeffs=[3.49548,0.0155532,-7.88013e-06,1.54339e-09,-1.07576e-13,-8239.2,13.6], Tmin=(832.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.0739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sH)"""),
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
    label = 'C#CC(F)[CH]F(5601)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u1 p0 c0 {2,S} {3,S} {8,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-15.2121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([323,546,651,945,1004,1329,3278,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.824138,'amu*angstrom^2'), symmetry=1, barrier=(18.9486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67582,'amu*angstrom^2'), symmetry=1, barrier=(38.5305,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77507,0.0167157,0.000142358,-4.50198e-07,3.93805e-10,-1827.44,11.3771], Tmin=(10,'K'), Tmax=(413.59,'K')), NASAPolynomial(coeffs=[6.12372,0.0275942,-1.89312e-05,6.1701e-09,-7.61207e-13,-2309.03,-1.35356], Tmin=(413.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-15.2121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C#CC(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C(CF)=C(F)CF(11293)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  C u2 p0 c0 {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-247.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932357,0.0733008,-8.16852e-05,5.69171e-08,-1.7355e-11,-29655.9,25.6055], Tmin=(100,'K'), Tmax=(776.338,'K')), NASAPolynomial(coeffs=[7.13526,0.0413388,-1.99256e-05,3.87841e-09,-2.7405e-13,-30619,-2.74911], Tmin=(776.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsCs) + group(CdCsCdF) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(CF)=C(F)[CH]F(9334)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {8,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u1 p0 c0 {3,S} {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-322.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916675,0.0714746,-7.6149e-05,4.36741e-08,-1.03283e-11,-38660.1,24.2303], Tmin=(100,'K'), Tmax=(1008.76,'K')), NASAPolynomial(coeffs=[11.4045,0.0298877,-1.43101e-05,2.80605e-09,-1.99939e-13,-40776,-26.4586], Tmin=(1008.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-322.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsCs) + group(CdCsCdF) + radical(Allyl_P) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = '[CH2]C(=CF)C(F)[CH]F(8562)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {8,D}
6  C u1 p0 c0 {2,S} {4,S} {10,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {5,D} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-276.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,350,440,435,1725,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.0102984,'amu*angstrom^2'), symmetry=1, barrier=(8.92065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.882121,'amu*angstrom^2'), symmetry=1, barrier=(20.2817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4557,'amu*angstrom^2'), symmetry=1, barrier=(56.4613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3166.93,'J/mol'), sigma=(5.33096,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=494.67 K, Pc=47.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668016,0.0775031,-9.26177e-05,5.93952e-08,-1.5481e-11,-33100.1,25.8646], Tmin=(100,'K'), Tmax=(927.264,'K')), NASAPolynomial(coeffs=[11.9727,0.0287377,-1.37324e-05,2.6803e-09,-1.90159e-13,-35196.6,-27.8204], Tmin=(927.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Csj(Cs-F1sCdH)(F1s)(H)) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(=CF)C(F)CF(9335)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {6,D} {12,S}
8  C u2 p0 c0 {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-250.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,194,682,905,1196,1383,3221,297.255,300.857,306.11,306.297,308.173],'cm^-1')),
        HinderedRotor(inertia=(0.768692,'amu*angstrom^2'), symmetry=1, barrier=(51.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792702,'amu*angstrom^2'), symmetry=1, barrier=(51.1049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786296,'amu*angstrom^2'), symmetry=1, barrier=(51.0542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972646,0.0716969,-6.95747e-05,3.87229e-08,-9.30951e-12,-30080.8,24.6734], Tmin=(100,'K'), Tmax=(969.68,'K')), NASAPolynomial(coeffs=[8.97806,0.0386746,-1.84935e-05,3.60475e-09,-2.55627e-13,-31633.4,-13.7018], Tmin=(969.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=CC(F)F)CF(11283)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {8,S}
7  C u0 p0 c0 {5,S} {6,D} {12,S}
8  C u2 p0 c0 {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-267.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12927,0.0692944,-7.05169e-05,4.50706e-08,-1.2958e-11,-32103.5,25.5766], Tmin=(100,'K'), Tmax=(808.791,'K')), NASAPolynomial(coeffs=[6.5494,0.0424885,-2.08021e-05,4.09211e-09,-2.91415e-13,-32980.3,0.577781], Tmin=(808.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'F[CH]C=C([CH]F)CF(9329)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {5,D} {8,S} {11,S}
7  C u1 p0 c0 {2,S} {5,S} {12,S}
8  C u1 p0 c0 {3,S} {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-318.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32791,0.0646391,-6.20157e-05,3.28911e-08,-7.51238e-12,-38206.9,24.9046], Tmin=(100,'K'), Tmax=(1012.92,'K')), NASAPolynomial(coeffs=[8.97403,0.0344452,-1.73032e-05,3.46331e-09,-2.49384e-13,-39755.9,-12.0818], Tmin=(1012.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Csj(Cd-CdH)(F1s)(H))"""),
)

species(
    label = '[CH]C(=C)C(F)C(F)F(9330)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {8,S}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  C u2 p0 c0 {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-275.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,2950,3100,1380,975,1025,1650,300.462,300.818,300.998,301.013,301.036],'cm^-1')),
        HinderedRotor(inertia=(0.790555,'amu*angstrom^2'), symmetry=1, barrier=(50.8417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.790252,'amu*angstrom^2'), symmetry=1, barrier=(50.8429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.790194,'amu*angstrom^2'), symmetry=1, barrier=(50.8408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.946238,0.0717905,-6.8573e-05,3.68033e-08,-8.45931e-12,-33025.6,24.7952], Tmin=(100,'K'), Tmax=(1014.44,'K')), NASAPolynomial(coeffs=[9.64354,0.0374963,-1.78636e-05,3.47805e-09,-2.46533e-13,-34790.2,-17.2889], Tmin=(1014.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-110.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (315.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (482.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (233.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-102.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-47.5807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (129.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (57.1721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (31.0789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1.44721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (14.9594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (294.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (157.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (25.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (81.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (34.2042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-41.3901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (23.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (46.2088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (78.7764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['CHFCHF[Z](58)', 'C#CCF(4379)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[CH]=CC(F)[CH]F(5593)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.44768e+28,'m^3/(mol*s)'), n=-6.4458, Ea=(84.8761,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.6591608693425184, var=5.4995193120720405, Tref=1000.0, N=19, data_mean=0.0, correlation='CH',), comment="""Estimated from node CH"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', '[CH]=C(F)C(F)[CH]F(5594)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(150.497,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]F(837)', '[CH]C(=CF)CF(11163)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['FCC1=CC(F)C1F(11304)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_1H;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['C=C(CF)C(F)=CF(9332)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', '[CH]=C(C=CF)CF(11142)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(67.3909,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CHFCHF[Z](58)', '[CH]=[C]CF(4380)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(2.01726,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH]=C(CF)C(F)=CF(11305)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH][CH]F(4261)', 'C#CCF(4379)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(11.0912,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2F(46)', 'C#CC(F)[CH]F(5601)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(19.4304,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[CH][CH]F(4261)', '[CH]=[C]CF(4380)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -18.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CHF(40)', '[CH]C(=CF)CF(11163)'],
    products = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['[CH]C(CF)=C(F)CF(11293)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.30951e+09,'s^-1'), n=1.14834, Ea=(136.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['[CH2]C(CF)=C(F)[CH]F(9334)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['[CH2]C(=CF)C(F)[CH]F(8562)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['[CH]C(=CF)C(F)CF(9335)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.67145,'s^-1'), n=3.19634, Ea=(69.5908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_1H;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['[CH]C(=CC(F)F)CF(11283)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(134.337,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['F[CH]C=C([CH]F)CF(9329)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(157.19,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(CF)C(F)[CH]F(9331)'],
    products = ['[CH]C(=C)C(F)C(F)F(9330)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(189.757,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #4027',
    isomers = [
        '[CH]=C(CF)C(F)[CH]F(9331)',
    ],
    reactants = [
        ('CHFCHF[Z](58)', 'C#CCF(4379)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4027',
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

