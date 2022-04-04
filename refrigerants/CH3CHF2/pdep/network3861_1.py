species(
    label = 'C=C(C)C(=O)OF(9774)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u0 p0 c0 {2,S} {3,D} {5,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-162.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,255.065,255.078,255.099,255.183,1470.18,4000],'cm^-1')),
        HinderedRotor(inertia=(0.138763,'amu*angstrom^2'), symmetry=1, barrier=(6.40797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138795,'amu*angstrom^2'), symmetry=1, barrier=(6.40793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671982,'amu*angstrom^2'), symmetry=1, barrier=(31.0315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41783,0.0633541,-9.21233e-05,8.43014e-08,-3.17228e-11,-19479.3,22.7631], Tmin=(100,'K'), Tmax=(775.048,'K')), NASAPolynomial(coeffs=[4.9767,0.0349231,-1.76221e-05,3.46501e-09,-2.4419e-13,-19728.7,8.45048], Tmin=(775.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
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
    label = 'C=C=C(6730)',
    structure = adjacencyList("""1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.00704622,2.78305e-05,-3.99443e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.706,'K')), NASAPolynomial(coeffs=[6.79957,0.00959978,-3.02068e-06,5.37826e-10,-3.92605e-14,19772.3,-12.7583], Tmin=(949.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=C(C)OF(580)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-54.0769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.164664,'amu*angstrom^2'), symmetry=1, barrier=(3.78595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0340763,'amu*angstrom^2'), symmetry=1, barrier=(7.51179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0696,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74756,0.0275218,-9.67407e-06,-2.80673e-09,1.90981e-12,-6504.87,9.64687], Tmin=(10,'K'), Tmax=(1181.53,'K')), NASAPolynomial(coeffs=[7.44792,0.0216801,-1.07455e-05,2.58688e-09,-2.44738e-13,-7845.96,-10.7974], Tmin=(1181.53,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-54.0769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDC(C)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=CC(=O)OF(1439)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 C u0 p0 c0 {4,D} {8,S} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-153.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.354365,'amu*angstrom^2'), symmetry=1, barrier=(8.14755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356401,'amu*angstrom^2'), symmetry=1, barrier=(8.19437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01811,0.0266015,-1.40043e-05,-2.86574e-08,3.81729e-11,-18381.3,13.9493], Tmin=(100,'K'), Tmax=(442.991,'K')), NASAPolynomial(coeffs=[4.78214,0.0179862,-9.59489e-06,1.97251e-09,-1.43646e-13,-18609.4,6.0652], Tmin=(442.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC(=C=O)COF(11537)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-131.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.423158,'amu*angstrom^2'), symmetry=1, barrier=(9.72923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425102,'amu*angstrom^2'), symmetry=1, barrier=(9.77393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423545,'amu*angstrom^2'), symmetry=1, barrier=(9.73813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885527,0.0760744,-0.000120995,1.06509e-07,-3.66089e-11,-15694.4,22.3521], Tmin=(100,'K'), Tmax=(851.592,'K')), NASAPolynomial(coeffs=[7.5438,0.0301212,-1.41974e-05,2.66236e-09,-1.80782e-13,-16296.1,-5.57551], Tmin=(851.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'C=C=C(OC)OF(11538)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {5,S}
4  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-6.71109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,348.566,349.103,349.482],'cm^-1')),
        HinderedRotor(inertia=(0.651957,'amu*angstrom^2'), symmetry=1, barrier=(14.9901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651796,'amu*angstrom^2'), symmetry=1, barrier=(14.9861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368076,'amu*angstrom^2'), symmetry=1, barrier=(14.9877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969516,0.0704851,-8.92078e-05,6.08747e-08,-1.67786e-11,-701.366,20.8563], Tmin=(100,'K'), Tmax=(881.851,'K')), NASAPolynomial(coeffs=[11.0216,0.0248903,-1.16535e-05,2.24555e-09,-1.5782e-13,-2474.28,-26.3755], Tmin=(881.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.71109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C([CH2])C(=O)OF(9827)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {3,D} {4,S}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (56.9789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.727545,0.0781347,-0.000116595,9.80255e-08,-3.33583e-11,6964.93,24.7336], Tmin=(100,'K'), Tmax=(779.347,'K')), NASAPolynomial(coeffs=[9.05164,0.0303983,-1.50685e-05,2.92494e-09,-2.04276e-13,5819.7,-12.3735], Tmin=(779.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.9789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJC(C)C=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(=C)C([O])OF(11539)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u1 p0 c0 {5,S} {9,S} {12,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (118.722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,627.506,627.695],'cm^-1')),
        HinderedRotor(inertia=(0.108083,'amu*angstrom^2'), symmetry=1, barrier=(2.48504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108005,'amu*angstrom^2'), symmetry=1, barrier=(2.48325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01654,'amu*angstrom^2'), symmetry=1, barrier=(46.3642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21633,0.0668132,-8.77724e-05,6.92041e-08,-2.30379e-11,14374,24.4059], Tmin=(100,'K'), Tmax=(723.561,'K')), NASAPolynomial(coeffs=[7.46895,0.0322473,-1.61145e-05,3.18053e-09,-2.25851e-13,13469.2,-3.73626], Tmin=(723.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(C)C([O])OF(11540)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
5  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u1 p0 c0 {6,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (214.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.134613,'amu*angstrom^2'), symmetry=1, barrier=(3.09502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134109,'amu*angstrom^2'), symmetry=1, barrier=(3.08343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.435654,'amu*angstrom^2'), symmetry=1, barrier=(10.0165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987103,0.0748657,-0.000121109,1.12052e-07,-4.08758e-11,25876.8,26.0318], Tmin=(100,'K'), Tmax=(815.916,'K')), NASAPolynomial(coeffs=[5.9675,0.0344801,-1.75047e-05,3.41101e-09,-2.37626e-13,25595.7,6.27491], Tmin=(815.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH2])=C(O)OF(11541)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {5,S} {12,S}
4  C u0 p0 c0 {5,D} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {3,S} {4,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {8,S} {11,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (61.9438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,757.103],'cm^-1')),
        HinderedRotor(inertia=(1.1698,'amu*angstrom^2'), symmetry=1, barrier=(26.896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.576066,'amu*angstrom^2'), symmetry=1, barrier=(16.1766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191069,'amu*angstrom^2'), symmetry=1, barrier=(78.0237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701174,'amu*angstrom^2'), symmetry=1, barrier=(16.1626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441984,0.0761526,-9.15895e-05,5.48732e-08,-1.27847e-11,7580.12,23.369], Tmin=(100,'K'), Tmax=(1056.87,'K')), NASAPolynomial(coeffs=[16.2926,0.016161,-6.44325e-06,1.16287e-09,-7.95309e-14,4229.75,-53.9778], Tmin=(1056.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(C)=C(O)OF(11542)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,S} {12,S}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {3,S} {5,D}
7  C u2 p0 c0 {5,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (129.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,316.23,316.241,316.242,316.245],'cm^-1')),
        HinderedRotor(inertia=(0.73298,'amu*angstrom^2'), symmetry=1, barrier=(52.0168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732993,'amu*angstrom^2'), symmetry=1, barrier=(52.0167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.733,'amu*angstrom^2'), symmetry=1, barrier=(52.0167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732984,'amu*angstrom^2'), symmetry=1, barrier=(52.0169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.54489,0.0797801,-0.000103513,7.45134e-08,-2.17111e-11,15711.8,23.5534], Tmin=(100,'K'), Tmax=(836.51,'K')), NASAPolynomial(coeffs=[11.0951,0.029332,-1.30521e-05,2.42001e-09,-1.65421e-13,13946.8,-25.4618], Tmin=(836.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(AllylJ2_triplet)"""),
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
    label = 'C=C(C)C([O])=O(11543)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {1,S} {2,D} {4,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-71.4828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180,2045.29,2045.53,2045.53],'cm^-1')),
        HinderedRotor(inertia=(0.143086,'amu*angstrom^2'), symmetry=1, barrier=(3.28983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143078,'amu*angstrom^2'), symmetry=1, barrier=(3.28965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0812,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21593,0.0487543,-7.96193e-05,8.70758e-08,-3.62256e-11,-8542.46,19.8924], Tmin=(100,'K'), Tmax=(820.266,'K')), NASAPolynomial(coeffs=[-0.730019,0.0396879,-2.01897e-05,3.94865e-09,-2.75876e-13,-7270.87,38.3264], Tmin=(820.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.4828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ)"""),
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
    label = 'C=C(C)[C]=O(9768)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  C u1 p0 c0 {1,D} {3,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (38.9217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.63819,'amu*angstrom^2'), symmetry=1, barrier=(14.6733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.637697,'amu*angstrom^2'), symmetry=1, barrier=(14.6619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0818,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84854,0.0323949,-1.85545e-05,4.53908e-09,-4.26922e-13,4714.82,14.4489], Tmin=(100,'K'), Tmax=(2340.37,'K')), NASAPolynomial(coeffs=[12.7643,0.0154475,-7.69257e-06,1.445e-09,-9.64098e-14,73.4875,-41.8206], Tmin=(2340.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.9217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'CH3(19)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=[C]C(=O)OF(2714)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 C u0 p0 c0 {6,D} {7,S} {8,S}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (90.8794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2950,3100,1380,975,1025,1650,1685,370,180,180,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.131767,'amu*angstrom^2'), symmetry=1, barrier=(3.02959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141118,'amu*angstrom^2'), symmetry=1, barrier=(3.24458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0451,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48131,0.0378621,-6.77965e-05,6.50273e-08,-2.46454e-11,10980.6,16.4924], Tmin=(100,'K'), Tmax=(754.612,'K')), NASAPolynomial(coeffs=[5.82747,0.014023,-8.2804e-06,1.73174e-09,-1.25721e-13,10649.3,2.44232], Tmin=(754.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.8794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=CJC=O)"""),
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
    label = '[CH2]C(=C)C(=O)OF(8746)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {5,S} {6,S} {7,D}
5  C u0 p0 c0 {2,S} {3,D} {4,S}
6  C u1 p0 c0 {4,S} {8,S} {9,S}
7  C u0 p0 c0 {4,D} {10,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-6.32872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,343.29,343.585,343.762,345.153,1885.27,4000],'cm^-1')),
        HinderedRotor(inertia=(0.36421,'amu*angstrom^2'), symmetry=1, barrier=(31.8866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0774783,'amu*angstrom^2'), symmetry=1, barrier=(6.46498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367253,'amu*angstrom^2'), symmetry=1, barrier=(31.9113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52787,0.0593847,-7.44525e-05,5.40583e-08,-1.64782e-11,-676.654,21.1921], Tmin=(100,'K'), Tmax=(787.99,'K')), NASAPolynomial(coeffs=[7.80506,0.0275195,-1.3793e-05,2.73687e-09,-1.95376e-13,-1665.9,-7.59586], Tmin=(787.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.32872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C=O)CJ)"""),
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
    label = 'C=[C]C(1569)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,D} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (239.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.586526,'amu*angstrom^2'), symmetry=1, barrier=(13.4854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2624,0.0120171,1.39877e-05,-2.15736e-08,7.7829e-12,28853.5,10.3012], Tmin=(100,'K'), Tmax=(1030.45,'K')), NASAPolynomial(coeffs=[5.02497,0.0154866,-6.07289e-06,1.11597e-09,-7.791e-14,27942.8,-0.91149], Tmin=(1030.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""propen2yl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C(C)C(=O)OF(8728)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u0 p0 c0 {2,S} {3,D} {5,S}
7  C u1 p0 c0 {5,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (84.4159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,229.645,232.276,232.44,242.537,1958.18],'cm^-1')),
        HinderedRotor(inertia=(0.00320431,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176998,'amu*angstrom^2'), symmetry=1, barrier=(6.40215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172964,'amu*angstrom^2'), symmetry=1, barrier=(6.41544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30499,0.0673385,-0.000110345,1.03492e-07,-3.82243e-11,10242.1,23.5707], Tmin=(100,'K'), Tmax=(811.346,'K')), NASAPolynomial(coeffs=[5.46781,0.0316504,-1.63288e-05,3.20398e-09,-2.24095e-13,10065.8,7.43375], Tmin=(811.346,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.4159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(C)C(=O)OF-2(11544)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {3,D} {4,S}
7  C u0 p1 c0 {4,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (85.5667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7731,0.052636,-4.14662e-05,1.69182e-08,-2.94364e-12,10368.1,21.8463], Tmin=(100,'K'), Tmax=(1286.89,'K')), NASAPolynomial(coeffs=[9.30488,0.0292251,-1.41783e-05,2.78185e-09,-1.97393e-13,8429.61,-16.3899], Tmin=(1286.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.5667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-188.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (23.6055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (221.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (34.4155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (261.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (31.4936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (133.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (229.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (48.1152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (115.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-46.9379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (93.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (178.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (157.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (155.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (247.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (64.1232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(C)C(=O)OF(9774)'],
    products = ['HF(38)', 'CO2(14)', 'C=C=C(6730)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00856353,'s^-1'), n=4.62568, Ea=(22.7801,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'C=C(C)OF(580)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(244.77,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(25)', 'C=CC(=O)OF(1439)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.11579e+08,'m^3/(mol*s)'), n=-0.607843, Ea=(3.45865,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CC(=C=O)COF(11537)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(214.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C=C(OC)OF(11538)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(316.063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_N-2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_N-2R!H->C"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C([CH2])C(=O)OF(9827)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=C)C([O])OF(11539)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C)C([O])OF(11540)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])=C(O)OF(11541)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]C(C)=C(O)OF(11542)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C=C(C)C([O])=O(11543)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.11128e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]F(126)', 'C=C(C)[C]=O(9768)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(19)', 'C=[C]C(=O)OF(2714)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', '[CH2]C(=C)C(=O)OF(8746)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.5958e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_Sp-4C=3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_Sp-4C=3C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]OF(167)', 'C=[C]C(1569)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[CH]=C(C)C(=O)OF(8728)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -5.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(C)C(=O)OF-2(11544)'],
    products = ['C=C(C)C(=O)OF(9774)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(26.9032,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

network(
    label = 'PDepNetwork #3861',
    isomers = [
        'C=C(C)C(=O)OF(9774)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'C=C=C(6730)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3861',
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

