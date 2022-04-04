species(
    label = 'O=C(O)C(F)(F)C(F)[CH]F(11992)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 C u1 p0 c0 {4,S} {8,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1068.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,348,427,405,1245,1236,1280,259,529,569,1128,1321,1390,3140,334,575,1197,1424,3202,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0515875,0.0968894,-0.000141608,1.10646e-07,-3.48955e-11,-128358,30.671], Tmin=(100,'K'), Tmax=(773.46,'K')), NASAPolynomial(coeffs=[12.3469,0.0327671,-1.72476e-05,3.45197e-09,-2.46247e-13,-130276,-25.9589], Tmin=(773.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1068.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    label = 'F[CH]C=C(F)F(1923)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u1 p0 c0 {3,S} {4,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-421.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.0575606,'amu*angstrom^2'), symmetry=1, barrier=(22.7411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0431,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2837.71,'J/mol'), sigma=(4.65594,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=443.24 K, Pc=63.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81319,0.0145998,0.000100931,-3.06247e-07,2.60885e-10,-50732.9,11.339], Tmin=(10,'K'), Tmax=(410.915,'K')), NASAPolynomial(coeffs=[4.51267,0.0276468,-1.91779e-05,6.2113e-09,-7.5793e-13,-50958.1,6.54677], Tmin=(410.915,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-421.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""F[CH]CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.56331e-09,3.13594e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88611e-07,1.21037e-10,-7.84012e-15,-14180.9,6.71043], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OC(F)(F)C(F)[CH]F(5414)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {11,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {6,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-929.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09216,'amu*angstrom^2'), symmetry=1, barrier=(25.1109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09322,'amu*angstrom^2'), symmetry=1, barrier=(25.1352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08873,'amu*angstrom^2'), symmetry=1, barrier=(25.0319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63088,0.0261579,0.000198074,-5.99931e-07,4.87251e-10,-111838,13.1613], Tmin=(10,'K'), Tmax=(452.864,'K')), NASAPolynomial(coeffs=[8.59939,0.0386199,-2.98398e-05,1.03341e-08,-1.31531e-12,-112866,-13.2525], Tmin=(452.864,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-929.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)(F)C(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(O)[C]F(3707)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {5,S}
5 C u0 p1 c0 {1,S} {4,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-238.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,262,1290,469.04,469.047,469.052,469.06,469.072],'cm^-1')),
        HinderedRotor(inertia=(0.000766243,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000766279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93204,0.0255218,-2.44778e-05,1.25722e-08,-2.7382e-12,-28605.1,14.1557], Tmin=(100,'K'), Tmax=(1066.48,'K')), NASAPolynomial(coeffs=[6.35429,0.0126859,-6.42393e-06,1.28648e-09,-9.26033e-14,-29335,-2.57488], Tmin=(1066.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdCsOs) + group(CJ2_singlet-FCO)"""),
)

species(
    label = 'CHF2-CHF(67)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-481.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.191999,'amu*angstrom^2'), symmetry=1, barrier=(4.41444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0324,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.798,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03194,0.0198962,-1.06494e-05,1.96763e-09,3.06708e-14,-57935.2,10.3725], Tmin=(10,'K'), Tmax=(1633.48,'K')), NASAPolynomial(coeffs=[15.0609,-0.000860555,2.67159e-06,-1.12649e-09,1.457e-13,-62372.2,-50.8003], Tmin=(1633.48,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-481.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""F[CH]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(F)C(F)F(1796)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-708.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,235,523,627,1123,1142,1372,1406,3097,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.148471,'amu*angstrom^2'), symmetry=1, barrier=(3.41363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622339,'amu*angstrom^2'), symmetry=1, barrier=(14.3088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.049,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2779.34,'J/mol'), sigma=(4.72663,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=434.13 K, Pc=59.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59036,0.0442234,6.97126e-05,-6.87999e-07,1.17292e-09,-85159.3,12.4166], Tmin=(10,'K'), Tmax=(241.349,'K')), NASAPolynomial(coeffs=[5.75627,0.0338039,-2.38726e-05,7.89282e-09,-9.81921e-13,-85338.1,3.50893], Tmin=(241.349,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-708.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""F[CH]C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OF(140)',
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
    label = 'O=C=C(F)C(F)[CH]F(12415)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {5,S} {10,S}
8  C u0 p0 c0 {4,D} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-417.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,145,326,398,834,1303,334,575,1197,1424,3202,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.102482,'amu*angstrom^2'), symmetry=1, barrier=(2.35625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00475,'amu*angstrom^2'), symmetry=1, barrier=(23.1012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.989576,0.0721066,-0.000106407,7.96475e-08,-2.2557e-11,-50087.8,24.9473], Tmin=(100,'K'), Tmax=(671.979,'K')), NASAPolynomial(coeffs=[10.8093,0.0214675,-1.08111e-05,2.11049e-09,-1.47908e-13,-51584,-19.8362], Tmin=(671.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cd(Cdd-Od)CF) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sCdH)(F1s)(H))"""),
)

species(
    label = 'CF2CO(76)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
"""),
    E0 = (-323.721,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(77.9917,'amu')),
        NonlinearRotor(inertia=([46.9849,128.952,175.938],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([205.962,243.776,382.471,461.225,682.222,822,1318.38,1476.6,2266.44],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.0175,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89388,0.00774418,6.875e-05,-2.11624e-07,1.80328e-10,-38932.8,8.78086], Tmin=(10,'K'), Tmax=(423.673,'K')), NASAPolynomial(coeffs=[5.02729,0.0134621,-9.62363e-06,3.16978e-09,-3.91588e-13,-39176.2,2.54717], Tmin=(423.673,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-323.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ODCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)[CH]F(982)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {8,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-459.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.469256,'amu*angstrom^2'), symmetry=1, barrier=(10.7891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781609,'amu*angstrom^2'), symmetry=1, barrier=(17.9707,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82583,0.0126617,0.000111215,-3.29958e-07,2.73521e-10,-55256.6,11.0522], Tmin=(10,'K'), Tmax=(431.912,'K')), NASAPolynomial(coeffs=[5.43113,0.0236036,-1.64171e-05,5.39472e-09,-6.67747e-13,-55636,1.86846], Tmin=(431.912,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-459.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""OC(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(OF)=C(F)C(F)[CH]F(12416)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {5,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {7,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-576.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,174,267,591,721,1107,1278,1348,3273,323,467,575,827,1418,350,440,435,1725,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.309219,'amu*angstrom^2'), symmetry=1, barrier=(7.10956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310044,'amu*angstrom^2'), symmetry=1, barrier=(7.12852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309483,'amu*angstrom^2'), symmetry=1, barrier=(7.11562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76789,'amu*angstrom^2'), symmetry=1, barrier=(40.6474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.875041,0.120125,-0.000215877,1.93186e-07,-6.60154e-11,-69154,32.6556], Tmin=(100,'K'), Tmax=(851.128,'K')), NASAPolynomial(coeffs=[12.6383,0.0319294,-1.69348e-05,3.28079e-09,-2.2491e-13,-70560.1,-25.1071], Tmin=(851.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-576.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsCsCs) + radical(Csj(Cs-F1sCdH)(F1s)(H))"""),
)

species(
    label = 'OC(OC(F)[CH]F)=C(F)F(12417)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u1 p0 c0 {2,S} {7,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-901.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,487,638,688,1119,1325,1387,3149,350,440,435,1725,334,575,1197,1424,3202,182,240,577,636,1210,1413,192.753,192.753,192.753,192.754],'cm^-1')),
        HinderedRotor(inertia=(0.836832,'amu*angstrom^2'), symmetry=1, barrier=(22.0631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836829,'amu*angstrom^2'), symmetry=1, barrier=(22.0632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836829,'amu*angstrom^2'), symmetry=1, barrier=(22.0632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836833,'amu*angstrom^2'), symmetry=1, barrier=(22.0632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.696944,0.107975,-0.000153093,1.05865e-07,-2.85257e-11,-108233,29.7135], Tmin=(100,'K'), Tmax=(913.545,'K')), NASAPolynomial(coeffs=[19.0604,0.0214654,-1.10469e-05,2.20494e-09,-1.57803e-13,-111843,-63.8175], Tmin=(913.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-901.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsCs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = '[CH]C(F)C(F)(F)C(=O)O(12418)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,S} {12,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u2 p0 c0 {7,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-638.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,348,427,405,1245,1236,1280,250,417,511,1155,1315,1456,3119,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59335,0.0808739,-0.000104748,7.18261e-08,-1.99758e-11,-76647.6,27.5308], Tmin=(100,'K'), Tmax=(871.128,'K')), NASAPolynomial(coeffs=[11.9752,0.028612,-1.47599e-05,2.95977e-09,-2.12612e-13,-78630.6,-25.8101], Tmin=(871.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-638.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(223)',
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
    label = 'O=C(O)C(F)(F)[CH]F(4018)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u1 p0 c0 {3,S} {6,S} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-837.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,249,341,381,651,928,1217,1251,334,575,1197,1424,3202,311.257,311.322,311.327,311.375,864.915,2604.63],'cm^-1')),
        HinderedRotor(inertia=(0.00871446,'amu*angstrom^2'), symmetry=1, barrier=(41.9529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609531,'amu*angstrom^2'), symmetry=1, barrier=(41.9528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609712,'amu*angstrom^2'), symmetry=1, barrier=(41.9535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2393,0.0659897,-9.09292e-05,6.82827e-08,-2.09053e-11,-100626,24.0799], Tmin=(100,'K'), Tmax=(793.192,'K')), NASAPolynomial(coeffs=[9.44075,0.0246279,-1.27054e-05,2.53267e-09,-1.80721e-13,-101927,-13.5867], Tmin=(793.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-837.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + group(Cds-OdCsOs) + radical(Csj(Cs-F1sF1sCO)(F1s)(H))"""),
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
    label = 'O=C(O)C(F)(F)C(F)[C]F(12419)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 C u2 p0 c0 {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-890.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,348,427,405,1245,1236,1280,250,417,511,1155,1315,1456,3119,90,150,250,247,323,377,431,433,1065,1274],'cm^-1')),
        HinderedRotor(inertia=(0.000942177,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29458,'amu*angstrom^2'), symmetry=1, barrier=(37.0732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160729,'amu*angstrom^2'), symmetry=1, barrier=(37.0557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230202,0.0904167,-0.00012957,9.80867e-08,-2.99777e-11,-106921,28.7538], Tmin=(100,'K'), Tmax=(796.834,'K')), NASAPolynomial(coeffs=[12.2202,0.0302283,-1.62679e-05,3.29262e-09,-2.36725e-13,-108832,-26.3679], Tmin=(796.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CsCCl_triplet)"""),
)

species(
    label = 'O[C]1OC(F)C(F)C1(F)F(12420)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {10,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1027.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.132312,0.0776822,-7.91024e-05,4.04631e-08,-7.91994e-12,-123375,27.549], Tmin=(100,'K'), Tmax=(1376.39,'K')), NASAPolynomial(coeffs=[19.184,0.0149912,-3.63776e-06,4.51033e-10,-2.3893e-14,-128072,-69.5563], Tmin=(1376.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1027.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(Tetrahydrofuran) + radical(Cs_P) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = '[O]C1(O)C(F)C(F)C1(F)F(12249)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
10 C u0 p0 c0 {2,S} {8,S} {9,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-954.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.105524,0.0795055,-7.90503e-05,3.77713e-08,-6.98223e-12,-114674,29.5179], Tmin=(100,'K'), Tmax=(1326.26,'K')), NASAPolynomial(coeffs=[21.245,0.0151122,-6.22127e-06,1.16256e-09,-8.14538e-14,-120337,-79.5145], Tmin=(1326.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-954.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(CC(C)(O)OJ) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'O=C(O)C(F)(F)C=CF(12421)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {12,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,D} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u0 p0 c0 {3,S} {7,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-959.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,271.573,271.574,271.576,271.579,271.579,271.585],'cm^-1')),
        HinderedRotor(inertia=(0.715312,'amu*angstrom^2'), symmetry=1, barrier=(37.4381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715337,'amu*angstrom^2'), symmetry=1, barrier=(37.4383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715352,'amu*angstrom^2'), symmetry=1, barrier=(37.4382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968995,0.0728733,-9.07439e-05,6.29646e-08,-1.81651e-11,-115265,25.4749], Tmin=(100,'K'), Tmax=(832.714,'K')), NASAPolynomial(coeffs=[9.59074,0.0314573,-1.61382e-05,3.23433e-09,-2.32341e-13,-116700,-14.5416], Tmin=(832.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-959.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(CdCFH)"""),
)

species(
    label = 'O=C(O)[C](F)F(1710)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {3,S} {4,D} {6,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-628.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,179,346,818,1406,1524,180,180,589.114,744.244,744.252],'cm^-1')),
        HinderedRotor(inertia=(0.0823368,'amu*angstrom^2'), symmetry=1, barrier=(32.3641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136253,'amu*angstrom^2'), symmetry=1, barrier=(53.5561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78509,0.0221402,1.13574e-05,-4.46075e-08,2.47599e-11,-75649.1,11.1505], Tmin=(10,'K'), Tmax=(745.082,'K')), NASAPolynomial(coeffs=[6.90617,0.0200223,-1.3848e-05,4.31269e-09,-5.01646e-13,-76520.5,-5.71514], Tmin=(745.082,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-628.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(O)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHFCHF[Z](59)',
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
    label = 'O=C(O)C(F)(F)C(F)=CF(8839)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 C u0 p0 c0 {4,S} {8,D} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1126.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,194,682,905,1196,1383,3221,180,180,180,180,180,779.714],'cm^-1')),
        HinderedRotor(inertia=(1.43465,'amu*angstrom^2'), symmetry=1, barrier=(32.9853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83876,'amu*angstrom^2'), symmetry=1, barrier=(42.2768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84118,'amu*angstrom^2'), symmetry=1, barrier=(42.3324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.139541,0.0923039,-0.000144126,1.2149e-07,-4.11474e-11,-135315,28.4785], Tmin=(100,'K'), Tmax=(757.37,'K')), NASAPolynomial(coeffs=[11.1312,0.0309875,-1.62211e-05,3.21236e-09,-2.26659e-13,-136886,-20.8774], Tmin=(757.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1126.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O=C(O)[C](F)C(F)[CH]F(12422)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {12,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {6,S} {8,S}
8  C u0 p0 c0 {4,S} {5,D} {7,S}
9  C u1 p0 c0 {3,S} {6,S} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-701.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,234,347,1316,1464,334,575,1197,1424,3202,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647611,0.0799326,-0.000104037,7.2779e-08,-2.07663e-11,-84275.4,28.0793], Tmin=(100,'K'), Tmax=(848.157,'K')), NASAPolynomial(coeffs=[11.2792,0.029792,-1.53605e-05,3.07631e-09,-2.20697e-13,-86078.8,-21.4609], Tmin=(848.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-701.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CsCOCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'O=C(O)C(F)(F)[CH][CH]F(12423)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {12,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u1 p0 c0 {3,S} {7,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-690.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,249,341,381,651,928,1217,1251,3025,407.5,1350,352.5,334,575,1197,1424,3202,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691245,0.0794078,-0.000109557,8.30159e-08,-2.57494e-11,-82984.2,30.3819], Tmin=(100,'K'), Tmax=(781.968,'K')), NASAPolynomial(coeffs=[10.2332,0.0305974,-1.59261e-05,3.18992e-09,-2.28268e-13,-84476.5,-13.3055], Tmin=(781.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-690.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + group(Cds-OdCsOs) + radical(CCJCC=O) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=[C]C(F)(F)C(F)[CH]F(12424)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u1 p0 c0 {4,S} {6,S} {11,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-643.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,169,350,401,508,725,1243,1223,334,575,1197,1424,3202,1855,455,950,296.264,296.952],'cm^-1')),
        HinderedRotor(inertia=(0.108967,'amu*angstrom^2'), symmetry=1, barrier=(6.77268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109373,'amu*angstrom^2'), symmetry=1, barrier=(6.78091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529482,'amu*angstrom^2'), symmetry=1, barrier=(32.9698,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631016,0.0828245,-0.000121527,8.6179e-08,-2.045e-11,-77300.7,28.2638], Tmin=(100,'K'), Tmax=(609.422,'K')), NASAPolynomial(coeffs=[11.1087,0.0266637,-1.43341e-05,2.87149e-09,-2.04048e-13,-78811.9,-19.0172], Tmin=(609.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(COj(Cs-F1sF1sCs)(O2d))"""),
)

species(
    label = '[O]C(=O)C(F)(F)C(F)[CH]F(12425)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {7,S} {12,S}
10 C u0 p0 c0 {5,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-842.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,231,348,427,405,1245,1236,1280,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0858923,0.0971007,-0.000165846,1.53863e-07,-5.61e-11,-101220,30.3254], Tmin=(100,'K'), Tmax=(800.131,'K')), NASAPolynomial(coeffs=[8.14613,0.0382932,-2.08943e-05,4.17189e-09,-2.94204e-13,-101918,-3.06008], Tmin=(800.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-842.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]F(493)',
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
    label = 'O=[C]O(158)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-192.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'F[CH]C(F)[C](F)F(919)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 C u1 p0 c0 {4,S} {5,S} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-506.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.165164,'amu*angstrom^2'), symmetry=1, barrier=(3.79745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16527,'amu*angstrom^2'), symmetry=1, barrier=(3.79988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23848,0.0729037,-0.000142734,1.40486e-07,-5.15257e-11,-60798.1,23.0025], Tmin=(100,'K'), Tmax=(853.927,'K')), NASAPolynomial(coeffs=[5.03329,0.027613,-1.48439e-05,2.90822e-09,-2.00632e-13,-60443,11.1679], Tmin=(853.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'O=C(O)C(F)(F)[C](F)[CH]F(12426)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {7,S} {10,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 C u1 p0 c0 {4,S} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-877.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,249,341,381,651,928,1217,1251,212,367,445,1450,334,575,1197,1424,3202,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20805,0.103589,-0.000180625,1.63493e-07,-5.7498e-11,-105447,33.2247], Tmin=(100,'K'), Tmax=(819.032,'K')), NASAPolynomial(coeffs=[10.3637,0.0334134,-1.81392e-05,3.58938e-09,-2.50885e-13,-106557,-11.8704], Tmin=(819.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-877.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    label = 'O=C(O)[C](F)C=CF(12427)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {11,S}
4  O u0 p2 c0 {7,D}
5  C u1 p0 c0 {1,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {3,S} {4,D} {5,S}
8  C u0 p0 c0 {2,S} {6,D} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-538.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,329.649,329.713,329.719,329.726,1149.63,2425.45],'cm^-1')),
        HinderedRotor(inertia=(0.496684,'amu*angstrom^2'), symmetry=1, barrier=(38.314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.49661,'amu*angstrom^2'), symmetry=1, barrier=(38.3138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.496457,'amu*angstrom^2'), symmetry=1, barrier=(38.3137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44207,0.0618999,-7.45502e-05,5.15452e-08,-1.50807e-11,-64668.9,23.9747], Tmin=(100,'K'), Tmax=(815.753,'K')), NASAPolynomial(coeffs=[8.02302,0.0296306,-1.52139e-05,3.05327e-09,-2.19599e-13,-65742.6,-6.43439], Tmin=(815.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(CdCFH) + radical(Cs_P)"""),
)

species(
    label = 'O=C(O)[C](F)C(F)=CF(12428)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u0 p2 c0 {8,D}
6  C u1 p0 c0 {2,S} {7,S} {8,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  C u0 p0 c0 {3,S} {7,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-711.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,271,519,563,612,1379,194,682,905,1196,1383,3221,224.03,224.04,224.046,969.696,1757.2,1757.2],'cm^-1')),
        HinderedRotor(inertia=(3.35876,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08754,'amu*angstrom^2'), symmetry=1, barrier=(38.74,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08768,'amu*angstrom^2'), symmetry=1, barrier=(38.7401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688934,0.0794299,-0.000122195,1.03847e-07,-3.57639e-11,-85480.9,27.0177], Tmin=(100,'K'), Tmax=(748.304,'K')), NASAPolynomial(coeffs=[9.48685,0.0292631,-1.53429e-05,3.04784e-09,-2.15703e-13,-86709.7,-12.2891], Tmin=(748.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cs_P)"""),
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
    label = 'O=C(O)C(F)(F)C(F)[C]F-2(12429)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 C u0 p1 c0 {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-895.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,348,427,405,1245,1236,1280,353,444,1253,3145,617,898,1187,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.426773,0.0860317,-0.000121047,9.16632e-08,-2.82657e-11,-107586,29.8258], Tmin=(100,'K'), Tmax=(787.446,'K')), NASAPolynomial(coeffs=[11.2019,0.031299,-1.67908e-05,3.40032e-09,-2.44726e-13,-109283,-19.5834], Tmin=(787.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-895.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCCFH) + group(Cds-OdCsOs) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'O=C(O)C(F)(F)[C](F)CF(12430)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,S} {13,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1070.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0420585,0.0950317,-0.000146027,1.24381e-07,-4.29919e-11,-128576,31.8055], Tmin=(100,'K'), Tmax=(745.897,'K')), NASAPolynomial(coeffs=[10.4411,0.0354748,-1.86353e-05,3.70744e-09,-2.62673e-13,-130022,-14.6081], Tmin=(745.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1070.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(=O)C(F)(F)C(F)CF(12431)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {10,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
10 C u0 p0 c0 {5,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-1034.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,231,348,427,405,1245,1236,1280,528,1116,1182,1331,1402,1494,3075,3110,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.992692,0.0781516,-7.71722e-05,2.78491e-09,3.82734e-11,-124378,26.6909], Tmin=(100,'K'), Tmax=(500.073,'K')), NASAPolynomial(coeffs=[8.03585,0.040724,-2.1626e-05,4.34975e-09,-3.11221e-13,-125318,-4.77132], Tmin=(500.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1034.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CCOJ)"""),
)

species(
    label = 'O=C(O)[C](F)C(F)C(F)F(12262)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {10,S} {13,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {10,S}
10 C u0 p0 c0 {5,S} {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1110.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.470285,0.0826836,-9.55007e-05,5.77099e-08,-1.42039e-11,-133404,29.0721], Tmin=(100,'K'), Tmax=(976.3,'K')), NASAPolynomial(coeffs=[13.1709,0.0306481,-1.55531e-05,3.118e-09,-2.24678e-13,-135884,-31.8966], Tmin=(976.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1110.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsOs) + radical(CsCOCsF1s)"""),
)

species(
    label = 'O=C(O)C(F)(F)[CH]C(F)F(12432)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {10,S} {13,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u0 p0 c0 {5,S} {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1104.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484199,0.084047,-0.000104392,6.92411e-08,-1.88125e-11,-132722,30.7786], Tmin=(100,'K'), Tmax=(886.516,'K')), NASAPolynomial(coeffs=[11.9514,0.0323058,-1.6844e-05,3.40345e-09,-2.45812e-13,-134756,-23.1625], Tmin=(886.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1104.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CCJCC=O)"""),
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
    E0 = (-410.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-347.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-63.9988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-314.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (128.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-62.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (13.1273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-340.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-93.8205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-150.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-206.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-543.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-483.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-350.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-455.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-442.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-157.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-146.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-143.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-159.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-228.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-226.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-194.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-68.6584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-330.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-227.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-208.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-460.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-536.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-390.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-432.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    products = ['HF(38)', 'CO2(14)', 'F[CH]C=C(F)F(1923)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.30828e+07,'s^-1'), n=1.83354, Ea=(185.991,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.49335582114639515, var=20.078174816455903, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_7R!H->O',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_7R!H->O"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'OC(F)(F)C(F)[CH]F(5414)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(229.712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C(O)[C]F(3707)', 'CHF2-CHF(67)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(184.241,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CO2(14)', 'F[CH]C(F)C(F)F(1796)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0347248,'m^3/(mol*s)'), n=2.50667, Ea=(324.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO2_Cdd;Cs_H] for rate rule [CO2_Cdd;C_ter]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_CO2"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(140)', 'O=C=C(F)C(F)[CH]F(12415)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000227174,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_Cdd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CF2CO(76)', 'OC(F)[CH]F(982)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.37684e-08,'m^3/(mol*s)'), n=3.46667, Ea=(249.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [doublebond;R_OH] + [Cd_Cdd;R_OR] for rate rule [Cd_Cdd;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC(OF)=C(F)C(F)[CH]F(12416)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(117.914,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OC(OC(F)[CH]F)=C(F)F(12417)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(89.2443,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[CH]C(F)C(F)(F)C(=O)O(12418)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]F(223)', 'O=C(O)C(F)(F)[CH]F(4018)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'O=C(O)C(F)(F)C(F)[C]F(12419)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    products = ['O[C]1OC(F)C(F)C1(F)F(12420)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.47116e+08,'s^-1'), n=0.669085, Ea=(53.4924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    products = ['[O]C1(O)C(F)C(F)C1(F)F(12249)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(113.622,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 110.7 to 113.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=C(O)C(F)(F)C=CF(12421)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(64.6448,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C(O)[C](F)F(1710)', 'CHFCHF[Z](59)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00504,'m^3/(mol*s)'), n=2.41, Ea=(12.2926,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H_Sp-2R!H=1R!H_Ext-4R!H-R_2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'O=C(O)C(F)(F)C(F)=CF(8839)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0579694,'m^3/(mol*s)'), n=2.57302, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.01164078162376979, var=0.9577230798162183, Tref=1000.0, N=13, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Sp-4CCl-1CCClOS"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'O=C(O)[C](F)C(F)[CH]F(12422)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -27.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', 'O=C(O)C(F)(F)[CH][CH]F(12423)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -29.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['OH(7)', 'O=[C]C(F)(F)C(F)[CH]F(12424)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(=O)C(F)(F)C(F)[CH]F(12425)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C(O)[C](F)F(1710)', 'F[CH][CH]F(493)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -8.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=[C]O(158)', 'F[CH]C(F)[C](F)F(919)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -8.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'O=C(O)C(F)(F)[C](F)[CH]F(12426)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0.261562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F2(78)', 'O=C(O)[C](F)C=CF(12427)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(7.00976,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'O=C(O)[C](F)C(F)=CF(12428)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(190.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CHF(40)', 'O=C(O)C(F)(F)[CH]F(4018)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(5)', 'O=C(O)C(F)(F)C(F)[C]F-2(12429)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(4.12257,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    products = ['O=C(O)C(F)(F)[C](F)CF(12430)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.30951e+09,'s^-1'), n=1.14834, Ea=(136.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(=O)C(F)(F)C(F)CF(12431)'],
    products = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.94594e+10,'s^-1'), n=0.343, Ea=(27.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_CCC;O_rad_out;Cs_H_out_1H] + [R5H_CCCC(O2d);O_rad_out;Cs_H_out] for rate rule [R5H_CCCC(O2d);O_rad_out;Cs_H_out_1H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    products = ['O=C(O)[C](F)C(F)C(F)F(12262)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(206.57,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=C(O)C(F)(F)C(F)[CH]F(11992)'],
    products = ['O=C(O)C(F)(F)[CH]C(F)F(12432)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(164.017,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3875',
    isomers = [
        'O=C(O)C(F)(F)C(F)[CH]F(11992)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'F[CH]C=C(F)F(1923)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3875',
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

