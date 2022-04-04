species(
    label = 'C=COOC=C(F)O[O](3537)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {3,S} {8,D} {11,S}
7  C u0 p0 c0 {2,S} {9,D} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {6,D}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-27.5597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.372343,'amu*angstrom^2'), symmetry=1, barrier=(8.56091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.371988,'amu*angstrom^2'), symmetry=1, barrier=(8.55275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372428,'amu*angstrom^2'), symmetry=1, barrier=(8.56285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17862,'amu*angstrom^2'), symmetry=1, barrier=(27.0988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416042,0.0859568,-0.000129021,1.09584e-07,-3.79308e-11,-3192.41,33.5712], Tmin=(100,'K'), Tmax=(753.286,'K')), NASAPolynomial(coeffs=[9.37781,0.0341569,-1.74849e-05,3.45013e-09,-2.43515e-13,-4423.05,-6.33195], Tmin=(753.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.5597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = '[O]OC(F)C=O(1537)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-328.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.18317,'amu*angstrom^2'), symmetry=1, barrier=(27.2033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603055,'amu*angstrom^2'), symmetry=1, barrier=(13.8654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3728.91,'J/mol'), sigma=(5.84942,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.45 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15439,0.0440801,-5.80197e-05,4.31461e-08,-1.32601e-11,-39422,18.3048], Tmin=(100,'K'), Tmax=(786.755,'K')), NASAPolynomial(coeffs=[7.18333,0.0185099,-9.26418e-06,1.82889e-09,-1.29953e-13,-40213.3,-4.7504], Tmin=(786.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'C=COC(F)(C=O)O[O](3593)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
7  C u0 p0 c0 {2,S} {9,D} {10,S}
8  C u0 p0 c0 {4,D} {6,S} {11,S}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-465.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.446579,0.0922083,-0.000110385,6.31645e-08,-1.38767e-11,-55844.5,29.2424], Tmin=(100,'K'), Tmax=(1125.47,'K')), NASAPolynomial(coeffs=[21.7043,0.0134827,-5.46174e-06,1.01356e-09,-7.11936e-14,-60830.6,-80.2408], Tmin=(1125.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-465.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=COCC=O(3538)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {5,S} {9,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {9,D} {13,S}
8  C u0 p0 c0 {4,D} {6,S} {12,S}
9  C u0 p0 c0 {1,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-342.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568171,0.079227,-9.32501e-05,5.69804e-08,-1.40286e-11,-41108.9,29.9021], Tmin=(100,'K'), Tmax=(982.044,'K')), NASAPolynomial(coeffs=[13.4427,0.0267876,-1.31536e-05,2.60685e-09,-1.86815e-13,-43637.6,-31.977], Tmin=(982.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-342.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = 'C=COO[CH]C(=O)F(3594)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u1 p0 c0 {2,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-261.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,611,648,830,1210,1753,2950,3100,1380,975,1025,1650,185.78,606.658],'cm^-1')),
        HinderedRotor(inertia=(0.926787,'amu*angstrom^2'), symmetry=1, barrier=(22.2421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.906128,'amu*angstrom^2'), symmetry=1, barrier=(22.216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.48201,'amu*angstrom^2'), symmetry=1, barrier=(58.1382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59903,'amu*angstrom^2'), symmetry=1, barrier=(38.3914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944161,0.0678435,-6.89303e-05,3.55729e-08,-7.40256e-12,-31386.1,27.2111], Tmin=(100,'K'), Tmax=(1150.29,'K')), NASAPolynomial(coeffs=[13.621,0.0237604,-1.1444e-05,2.25536e-09,-1.613e-13,-34302.4,-35.7222], Tmin=(1150.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(COCsFO) + group(Cds-CdsHH) + radical(OCJC=O)"""),
)

species(
    label = 'FC1=COOO1(3394)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u0 p0 c0 {2,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-63.0576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45648,-0.00296079,7.42328e-05,-9.36209e-08,3.49566e-11,-7551.03,17.4583], Tmin=(100,'K'), Tmax=(965.686,'K')), NASAPolynomial(coeffs=[9.89799,0.00855895,-2.99883e-06,6.60762e-10,-5.63624e-14,-10576.4,-22.616], Tmin=(965.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.0576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclopentane)"""),
)

species(
    label = 'CH2CHO(35)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,531.349,531.873,974.177,1641.05,1641.53],'cm^-1')),
        HinderedRotor(inertia=(0.00117252,'amu*angstrom^2'), symmetry=1, barrier=(2.241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=COOC1OO[C]1F(3595)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {2,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
7  C u1 p0 c0 {1,S} {4,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-67.3832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901081,0.070304,-6.85878e-05,3.4126e-08,-6.92923e-12,-7994.57,30.7073], Tmin=(100,'K'), Tmax=(1165.78,'K')), NASAPolynomial(coeffs=[13.336,0.0276368,-1.36874e-05,2.72999e-09,-1.96306e-13,-10893.8,-31.1915], Tmin=(1165.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.3832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(CsCFHO) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'FC1=COO[CH]COO1(3596)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u0 p2 c0 {3,S} {8,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {3,S} {6,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {1,S} {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-19.1546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64425,0.0467986,-1.38172e-05,-5.88432e-09,2.71338e-12,-2215.06,24.0779], Tmin=(100,'K'), Tmax=(1451.14,'K')), NASAPolynomial(coeffs=[13.1993,0.031856,-1.58491e-05,3.0785e-09,-2.14004e-13,-7348.96,-42.1053], Tmin=(1451.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.1546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFO) + ring(Cyclooctane) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'C=COO[CH]C1(F)OO1(3597)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u1 p0 c0 {4,S} {6,S} {10,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-108.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.588788,0.0770231,-8.15449e-05,4.42094e-08,-9.68397e-12,-12941.7,29.4011], Tmin=(100,'K'), Tmax=(1094.86,'K')), NASAPolynomial(coeffs=[14.2331,0.027174,-1.3249e-05,2.62307e-09,-1.88059e-13,-15929.4,-37.6611], Tmin=(1094.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(O2s-O2s-Cs(F)(C)) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]C1OOC=C(F)OO1(3598)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {3,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
7  C u1 p0 c0 {6,S} {12,S} {13,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {1,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-47.6667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89014,0.025432,0.000127702,-2.17755e-07,9.45998e-11,-5581.01,30.6665], Tmin=(100,'K'), Tmax=(919.244,'K')), NASAPolynomial(coeffs=[35.9676,-0.0156801,1.28052e-05,-2.44877e-09,1.50997e-13,-16741.9,-161.238], Tmin=(919.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.6667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + ring(Cycloheptane) + radical(CJCOOH) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    label = 'C=COOC=[C]O[O](3599)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {7,D} {9,S}
6  C u0 p0 c0 {2,S} {8,D} {12,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  C u1 p0 c0 {3,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (370.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.602258,'amu*angstrom^2'), symmetry=1, barrier=(13.8471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602241,'amu*angstrom^2'), symmetry=1, barrier=(13.8467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602848,'amu*angstrom^2'), symmetry=1, barrier=(13.8607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60252,'amu*angstrom^2'), symmetry=1, barrier=(13.8531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.776492,0.0721374,-8.83126e-05,5.56192e-08,-1.38541e-11,44687.6,33.3953], Tmin=(100,'K'), Tmax=(981.766,'K')), NASAPolynomial(coeffs=[13.5393,0.020138,-8.86501e-06,1.67052e-09,-1.16439e-13,42181.5,-27.9431], Tmin=(981.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]O[C](F)C=O(1584)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-189.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,280,501,1494,1531,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(4.16547,'amu*angstrom^2'), symmetry=1, barrier=(95.7722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.22148,'amu*angstrom^2'), symmetry=1, barrier=(97.0601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35167,0.0412352,-6.6895e-05,6.24069e-08,-2.28565e-11,-22685.8,17.8037], Tmin=(100,'K'), Tmax=(824.824,'K')), NASAPolynomial(coeffs=[4.73144,0.0198291,-1.00257e-05,1.94134e-09,-1.3458e-13,-22742.8,8.81534], Tmin=(824.824,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]OC=C(F)O[O](3398)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u1 p2 c0 {2,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {2,S} {7,D} {8,S}
7 C u0 p0 c0 {1,S} {3,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (19.1511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.233153,'amu*angstrom^2'), symmetry=1, barrier=(5.36065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234516,'amu*angstrom^2'), symmetry=1, barrier=(5.39199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68886,0.0630752,-0.000129038,1.29573e-07,-4.76594e-11,2374.86,24.991], Tmin=(100,'K'), Tmax=(871.403,'K')), NASAPolynomial(coeffs=[3.58325,0.025133,-1.33819e-05,2.57542e-09,-1.74843e-13,3155.11,22.4838], Tmin=(871.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.1511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = 'C2H3(29)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=C(F)O[O](1081)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u1 p0 c0 {4,D} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (157.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,293,496,537,1218,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.441569,'amu*angstrom^2'), symmetry=1, barrier=(10.1525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00596,0.051545,-0.000104345,9.71572e-08,-3.31603e-11,19053.1,15.8311], Tmin=(100,'K'), Tmax=(897.319,'K')), NASAPolynomial(coeffs=[7.11388,0.0105561,-5.37028e-06,9.95721e-10,-6.48097e-14,18869.9,-4.17116], Tmin=(897.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'C=CO[O](287)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (100.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.13511,'amu*angstrom^2'), symmetry=1, barrier=(26.0985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.83227,0.0191966,3.62976e-06,-2.04121e-08,9.62936e-12,12153.9,12.3877], Tmin=(100,'K'), Tmax=(970.344,'K')), NASAPolynomial(coeffs=[9.43949,0.00760096,-2.62382e-06,4.96063e-10,-3.7262e-14,10135.3,-23.0837], Tmin=(970.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.0012157,5.31615e-06,-4.8944e-09,1.45844e-12,-1038.59,4.68369], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15383,0.00167803,-7.69968e-07,1.51274e-10,-1.08781e-14,-1040.82,6.16752], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=COOC=[C]F(3600)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (157.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.718823,'amu*angstrom^2'), symmetry=1, barrier=(16.5271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.71891,'amu*angstrom^2'), symmetry=1, barrier=(16.5292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14813,'amu*angstrom^2'), symmetry=1, barrier=(26.3978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37223,0.0606733,-6.79754e-05,4.02956e-08,-9.72497e-12,18988.3,25.2957], Tmin=(100,'K'), Tmax=(996.072,'K')), NASAPolynomial(coeffs=[10.7932,0.0228408,-1.10028e-05,2.164e-09,-1.54465e-13,17111.5,-20.1179], Tmin=(996.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    label = 'C=COO[C]=C(F)O[O](3601)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,D} {10,S}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (212.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,293,496,537,1218,1685,370,180,1278.78],'cm^-1')),
        HinderedRotor(inertia=(0.297539,'amu*angstrom^2'), symmetry=1, barrier=(6.84102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298062,'amu*angstrom^2'), symmetry=1, barrier=(6.85303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299579,'amu*angstrom^2'), symmetry=1, barrier=(6.88791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300568,'amu*angstrom^2'), symmetry=1, barrier=(6.91065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.536935,0.088835,-0.000160122,1.54164e-07,-5.66946e-11,25632.4,36.2607], Tmin=(100,'K'), Tmax=(836.494,'K')), NASAPolynomial(coeffs=[5.6082,0.037622,-1.99376e-05,3.90643e-09,-2.71177e-13,25727.3,18.3388], Tmin=(836.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OOC=C(F)O[O](3602)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {6,D}
8  C u0 p0 c0 {9,D} {11,S} {12,S}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (212.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.274286,'amu*angstrom^2'), symmetry=1, barrier=(6.30637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274006,'amu*angstrom^2'), symmetry=1, barrier=(6.29994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274049,'amu*angstrom^2'), symmetry=1, barrier=(6.30092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27428,'amu*angstrom^2'), symmetry=1, barrier=(6.30623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.536935,0.088835,-0.000160122,1.54164e-07,-5.66946e-11,25632.4,36.2607], Tmin=(100,'K'), Tmax=(836.494,'K')), NASAPolynomial(coeffs=[5.6082,0.037622,-1.99376e-05,3.90643e-09,-2.71177e-13,25727.3,18.3388], Tmin=(836.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COOC=C(F)O[O](3603)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {6,D}
8  C u0 p0 c0 {3,S} {9,D} {11,S}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (219.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.439671,'amu*angstrom^2'), symmetry=1, barrier=(10.1089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439706,'amu*angstrom^2'), symmetry=1, barrier=(10.1097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08682,'amu*angstrom^2'), symmetry=1, barrier=(24.9882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439955,'amu*angstrom^2'), symmetry=1, barrier=(10.1154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277669,0.0902799,-0.000148621,1.30877e-07,-4.54874e-11,26530.1,34.468], Tmin=(100,'K'), Tmax=(803.368,'K')), NASAPolynomial(coeffs=[9.92653,0.0307786,-1.61273e-05,3.17335e-09,-2.22075e-13,25349.5,-7.6681], Tmin=(803.368,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[O]OC(F)=C=O(3604)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u1 p2 c0 {2,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (-105.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (91.0179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37136,0.043506,-8.64301e-05,8.52867e-08,-3.11222e-11,-12657.6,18.7731], Tmin=(100,'K'), Tmax=(866.808,'K')), NASAPolynomial(coeffs=[4.28985,0.0166729,-8.88159e-06,1.71369e-09,-1.1673e-13,-12314.7,13.688], Tmin=(866.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = 'CH3CHO(36)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.19,1305.32,1306.47,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136162,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'C=COOC#CO[O](3605)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,D} {9,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {8,T}
8  C u0 p0 c0 {3,S} {7,T}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (421.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2100,2250,500,550,180],'cm^-1')),
        HinderedRotor(inertia=(1.26988,'amu*angstrom^2'), symmetry=1, barrier=(29.197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26923,'amu*angstrom^2'), symmetry=1, barrier=(29.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27368,'amu*angstrom^2'), symmetry=1, barrier=(29.2844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2681,'amu*angstrom^2'), symmetry=1, barrier=(29.156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88,0.0717641,-8.89847e-05,5.63882e-08,-1.42504e-11,50838.2,29.1155], Tmin=(100,'K'), Tmax=(962.588,'K')), NASAPolynomial(coeffs=[12.9858,0.0214589,-1.05943e-05,2.09674e-09,-1.49999e-13,48507.6,-28.8264], Tmin=(962.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(O2s-OsCt) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'C=COO[C]=C(F)OO(3606)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u0 p2 c0 {3,S} {13,S}
6  C u0 p0 c0 {2,S} {8,D} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u0 p0 c0 {6,D} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (60.1798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,293,496,537,1218,2950,3100,1380,975,1025,1650,1685,370,180,1783.24],'cm^-1')),
        HinderedRotor(inertia=(0.273781,'amu*angstrom^2'), symmetry=1, barrier=(6.29476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273866,'amu*angstrom^2'), symmetry=1, barrier=(6.29671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71156,'amu*angstrom^2'), symmetry=1, barrier=(39.352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273999,'amu*angstrom^2'), symmetry=1, barrier=(6.29977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71116,'amu*angstrom^2'), symmetry=1, barrier=(39.3428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.295903,0.0921734,-0.000154684,1.44321e-07,-5.30068e-11,7360.95,36.3664], Tmin=(100,'K'), Tmax=(803.342,'K')), NASAPolynomial(coeffs=[7.0913,0.0394678,-2.10383e-05,4.17352e-09,-2.93521e-13,6878.04,8.86036], Tmin=(803.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.1798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OOC=C(F)OO(3607)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u0 p2 c0 {3,S} {13,S}
6  C u0 p0 c0 {2,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {6,D}
8  C u0 p0 c0 {9,D} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (60.1798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.22056,'amu*angstrom^2'), symmetry=1, barrier=(5.07112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57105,'amu*angstrom^2'), symmetry=1, barrier=(36.1215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22077,'amu*angstrom^2'), symmetry=1, barrier=(5.07593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220858,'amu*angstrom^2'), symmetry=1, barrier=(5.07795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.295903,0.0921734,-0.000154684,1.44321e-07,-5.30068e-11,7360.95,36.3664], Tmin=(100,'K'), Tmax=(803.342,'K')), NASAPolynomial(coeffs=[7.0913,0.0394678,-2.10383e-05,4.17352e-09,-2.93521e-13,6878.04,8.86036], Tmin=(803.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.1798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COOC=C(F)OO(3608)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {13,S}
6  C u0 p0 c0 {2,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {6,D}
8  C u0 p0 c0 {3,S} {9,D} {11,S}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (67.5317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.508395,'amu*angstrom^2'), symmetry=1, barrier=(11.689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507658,'amu*angstrom^2'), symmetry=1, barrier=(11.6721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.509111,'amu*angstrom^2'), symmetry=1, barrier=(11.7055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67839,'amu*angstrom^2'), symmetry=1, barrier=(38.5896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6754,'amu*angstrom^2'), symmetry=1, barrier=(38.5207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184146,0.0916513,-0.000135116,1.08615e-07,-3.54905e-11,8252.42,34.0585], Tmin=(100,'K'), Tmax=(745.543,'K')), NASAPolynomial(coeffs=[11.0691,0.0332487,-1.76076e-05,3.53359e-09,-2.52369e-13,6629.46,-15.2582], Tmin=(745.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.5317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=COOC=[C]OOF(3609)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {1,S} {4,S}
6  C u0 p0 c0 {2,S} {8,D} {10,S}
7  C u0 p0 c0 {3,S} {9,D} {11,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (313.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.06001,'amu*angstrom^2'), symmetry=1, barrier=(24.3717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06693,'amu*angstrom^2'), symmetry=1, barrier=(24.5308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05404,'amu*angstrom^2'), symmetry=1, barrier=(24.2345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267957,0.0839846,-0.000101146,6.13289e-08,-1.47167e-11,37885.4,35.8465], Tmin=(100,'K'), Tmax=(1015.79,'K')), NASAPolynomial(coeffs=[15.7721,0.0229321,-1.09908e-05,2.15948e-09,-1.54262e-13,34735.7,-39.1947], Tmin=(1015.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
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
    E0 = (-57.2492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-8.87781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (23.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-71.4055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-25.3198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (46.1053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-46.9296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-41.6374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (28.7757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (390.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-80.1302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (252.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (205.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (95.9134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (371.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (371.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (378.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-34.7146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (374.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (51.9178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (116.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (148.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (354.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['CH2CO(28)', '[O]OC(F)C=O(1537)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(22.881,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['C=COC(F)(C=O)O[O](3593)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(71.2524,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['[O]OC(F)=COCC=O(3538)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(104.014,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'C=COO[CH]C(=O)F(3594)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['FC1=COOO1(3394)', 'CH2CHO(35)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['C=COOC1OO[C]1F(3595)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['FC1=COO[CH]COO1(3596)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.56799e+08,'s^-1'), n=0.604362, Ea=(33.2006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['C=COO[CH]C1(F)OO1(3597)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(38.4928,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['[CH2]C1OOC=C(F)OO1(3598)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.19557e+09,'s^-1'), n=0.501475, Ea=(108.906,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'C=COOC=[C]O[O](3599)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]O[C](F)C=O(1584)', 'CH2CHO(35)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(160.285,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Ea raised from 157.2 to 160.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC=C(F)O[O](3398)', 'C2H3(29)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3e+06,'m^3/(mol*s)'), n=-6.55423e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_N-Sp-4R!H-2C_4R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_N-Sp-4R!H-2C_4R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(F)O[O](1081)', 'C=CO[O](287)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)', 'C=COOC=[C]F(3600)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'C=COO[C]=C(F)O[O](3601)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'C=[C]OOC=C(F)O[O](3602)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[CH]=COOC=C(F)O[O](3603)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=COOC=C(F)O[O](3537)'],
    products = ['[O]OC(F)=C=O(3604)', 'CH3CHO(36)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(45.4156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'C=COOC#CO[O](3605)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(286.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=COO[C]=C(F)OO(3606)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]OOC=C(F)OO(3607)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.81284e+06,'s^-1'), n=1.31318, Ea=(108.582,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R7H_OOCs4;Y_rad_out;XH_out] for rate rule [R7H_OOCs4;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.1622776601683795
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=COOC=C(F)OO(3608)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.21044e+06,'s^-1'), n=1.67676, Ea=(133.491,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R8H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=COOC=[C]OOF(3609)'],
    products = ['C=COOC=C(F)O[O](3537)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.8166,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1243',
    isomers = [
        'C=COOC=C(F)O[O](3537)',
    ],
    reactants = [
        ('CH2CO(28)', '[O]OC(F)C=O(1537)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1243',
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

