species(
    label = 'COOC=C(F)O[O](3351)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-132.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.313339,'amu*angstrom^2'), symmetry=1, barrier=(7.20429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31344,'amu*angstrom^2'), symmetry=1, barrier=(7.20659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.313344,'amu*angstrom^2'), symmetry=1, barrier=(7.2044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31344,'amu*angstrom^2'), symmetry=1, barrier=(7.2066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916009,0.0794777,-0.000138687,1.35832e-07,-5.12095e-11,-15888,29.2589], Tmin=(100,'K'), Tmax=(826.679,'K')), NASAPolynomial(coeffs=[4.00386,0.0395013,-2.07235e-05,4.0676e-09,-2.83567e-13,-15543,20.1236], Tmin=(826.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=COO(3391)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {2,S} {9,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {2,S} {7,D} {8,S}
7 C u0 p0 c0 {1,S} {3,S} {6,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-132.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.363982,'amu*angstrom^2'), symmetry=1, barrier=(8.36865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363921,'amu*angstrom^2'), symmetry=1, barrier=(8.36727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363788,'amu*angstrom^2'), symmetry=1, barrier=(8.3642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4277,0.0666494,-0.000124382,1.20591e-07,-4.42134e-11,-15895.7,25.1687], Tmin=(100,'K'), Tmax=(844.524,'K')), NASAPolynomial(coeffs=[5.19556,0.0267509,-1.43476e-05,2.81e-09,-1.94452e-13,-15745.7,12.2837], Tmin=(844.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = 'COC(F)(C=O)O[O](3392)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
7  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {4,D} {6,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-516.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.560131,0.0850201,-0.000142465,1.29475e-07,-4.58688e-11,-61946.5,26.0976], Tmin=(100,'K'), Tmax=(829.091,'K')), NASAPolynomial(coeffs=[7.79462,0.03283,-1.67667e-05,3.25275e-09,-2.25249e-13,-62552,-3.86556], Tmin=(829.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-516.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = 'COO[CH]C(=O)F(3393)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {7,S} {11,S}
7  C u0 p0 c0 {1,S} {4,D} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-367.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(0.0456809,'amu*angstrom^2'), symmetry=1, barrier=(46.6017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407826,'amu*angstrom^2'), symmetry=1, barrier=(9.37672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02764,'amu*angstrom^2'), symmetry=1, barrier=(46.6195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50731,0.0603782,-7.41736e-05,5.51351e-08,-1.7642e-11,-44083.9,22.6894], Tmin=(100,'K'), Tmax=(745.393,'K')), NASAPolynomial(coeffs=[6.88155,0.0315376,-1.61345e-05,3.22463e-09,-2.311e-13,-44885.1,-1.65887], Tmin=(745.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsHHH) + group(COCsFO) + radical(OCJC=O)"""),
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
    label = 'COOC1OO[C]1F(3395)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {3,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {5,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-172.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39137,0.0636717,-7.66777e-05,5.74355e-08,-1.8871e-11,-20689.2,26.4491], Tmin=(100,'K'), Tmax=(721.572,'K')), NASAPolynomial(coeffs=[6.43503,0.0357138,-1.85621e-05,3.74483e-09,-2.7002e-13,-21417.1,3.76196], Tmin=(721.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(Cs-OsHHH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'COO[CH]C1(F)OO1(3396)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
7  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {12,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-213.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08707,0.0703645,-8.97936e-05,6.79203e-08,-2.18065e-11,-25636.9,25.1091], Tmin=(100,'K'), Tmax=(746.909,'K')), NASAPolynomial(coeffs=[7.79584,0.0344371,-1.76429e-05,3.52223e-09,-2.52042e-13,-26639.1,-5.2993], Tmin=(746.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cs-OsHHH) + ring(O2s-O2s-Cs(F)(C)) + radical(CCsJOOC)"""),
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
    label = 'COOC=[C]O[O](3397)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {7,D} {11,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (265.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.391769,'amu*angstrom^2'), symmetry=1, barrier=(9.00753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524763,'amu*angstrom^2'), symmetry=1, barrier=(12.0653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.391003,'amu*angstrom^2'), symmetry=1, barrier=(8.98992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.391673,'amu*angstrom^2'), symmetry=1, barrier=(9.00534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23989,0.0661074,-9.97251e-05,8.46387e-08,-2.86907e-11,31993.6,29.2137], Tmin=(100,'K'), Tmax=(818.814,'K')), NASAPolynomial(coeffs=[7.94352,0.0258845,-1.23465e-05,2.34736e-09,-1.61545e-13,31146.4,-0.257306], Tmin=(818.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
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
    label = 'CO[O](183)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-0.245543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,891.071,1078.46],'cm^-1')),
        HinderedRotor(inertia=(0.29694,'amu*angstrom^2'), symmetry=1, barrier=(6.82724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.13,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.9717,-0.00529357,4.77334e-05,-5.77066e-08,2.2222e-11,-129.022,2.81501], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.5553,0.00912236,-3.23852e-06,5.18714e-10,-3.08834e-14,-1035.69,-3.99159], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-0.245543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3O2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'COOC=[C]F(3399)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u1 p0 c0 {1,S} {5,D}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (51.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(0.355257,'amu*angstrom^2'), symmetry=1, barrier=(8.16806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355106,'amu*angstrom^2'), symmetry=1, barrier=(8.16458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25302,'amu*angstrom^2'), symmetry=1, barrier=(51.8015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71556,0.0561569,-8.52205e-05,7.7858e-08,-2.87556e-11,6299.53,21.539], Tmin=(100,'K'), Tmax=(791.547,'K')), NASAPolynomial(coeffs=[5.2994,0.0284109,-1.43816e-05,2.81644e-09,-1.9754e-13,6034.03,6.99355], Tmin=(791.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    label = '[CH2]OOC=C(F)O[O](3400)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,D} {9,S}
7  C u0 p0 c0 {1,S} {4,S} {6,D}
8  C u1 p0 c0 {3,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (60.6674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.228653,'amu*angstrom^2'), symmetry=1, barrier=(5.25719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228129,'amu*angstrom^2'), symmetry=1, barrier=(5.24514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227627,'amu*angstrom^2'), symmetry=1, barrier=(5.23359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.3956,'amu*angstrom^2'), symmetry=1, barrier=(55.0795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863951,0.0797905,-0.00014055,1.34661e-07,-4.97505e-11,7399.04,29.1427], Tmin=(100,'K'), Tmax=(825.978,'K')), NASAPolynomial(coeffs=[5.5796,0.0347711,-1.8509e-05,3.64428e-09,-2.5423e-13,7376.74,11.8746], Tmin=(825.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.6674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(CsJOOC)"""),
)

species(
    label = 'COO[C]=C(F)O[O](3401)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {4,S} {8,D}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (106.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,293,496,537,1218,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.232379,'amu*angstrom^2'), symmetry=1, barrier=(5.34286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232408,'amu*angstrom^2'), symmetry=1, barrier=(5.34352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232593,'amu*angstrom^2'), symmetry=1, barrier=(5.34777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181053,'amu*angstrom^2'), symmetry=1, barrier=(5.36144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10712,0.0814608,-0.000166354,1.75572e-07,-6.77747e-11,12933.8,31.7011], Tmin=(100,'K'), Tmax=(853.751,'K')), NASAPolynomial(coeffs=[-0.00871653,0.0434031,-2.34383e-05,4.58759e-09,-3.16627e-13,14701.9,46.1467], Tmin=(853.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
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
    label = 'COOC#CO[O](3402)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {3,S} {6,T}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (316.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550],'cm^-1')),
        HinderedRotor(inertia=(0.863504,'amu*angstrom^2'), symmetry=1, barrier=(19.8537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91269,'amu*angstrom^2'), symmetry=1, barrier=(43.9765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863413,'amu*angstrom^2'), symmetry=1, barrier=(19.8516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9146,'amu*angstrom^2'), symmetry=1, barrier=(44.0205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23316,0.0671416,-0.00010589,9.35182e-08,-3.30816e-11,38148.9,25.3227], Tmin=(100,'K'), Tmax=(790.498,'K')), NASAPolynomial(coeffs=[7.54925,0.0269249,-1.39101e-05,2.73375e-09,-1.9176e-13,37408.3,-2.03217], Tmin=(790.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-OsHHH) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'COO[C]=C(F)OO(3403)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {3,S} {12,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-45.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,293,496,537,1218,1685,370,1367.65],'cm^-1')),
        HinderedRotor(inertia=(0.218075,'amu*angstrom^2'), symmetry=1, barrier=(5.01398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218061,'amu*angstrom^2'), symmetry=1, barrier=(5.01364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218087,'amu*angstrom^2'), symmetry=1, barrier=(5.01424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218292,'amu*angstrom^2'), symmetry=1, barrier=(5.01897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.29838,'amu*angstrom^2'), symmetry=1, barrier=(52.8442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84691,0.0850323,-0.000161744,1.66764e-07,-6.44772e-11,-5336.78,31.8749], Tmin=(100,'K'), Tmax=(835.117,'K')), NASAPolynomial(coeffs=[1.56911,0.0450808,-2.4439e-05,4.83054e-09,-3.36936e-13,-4184.87,36.1397], Tmin=(835.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]OOC=C(F)OO(3404)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {12,S}
6  C u0 p0 c0 {2,S} {7,D} {9,S}
7  C u0 p0 c0 {1,S} {4,S} {6,D}
8  C u1 p0 c0 {3,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-91.3373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633589,0.0829993,-0.000134655,1.24266e-07,-4.58728e-11,-10872.8,29.2105], Tmin=(100,'K'), Tmax=(782.543,'K')), NASAPolynomial(coeffs=[7.00686,0.0367169,-1.96695e-05,3.92588e-09,-2.77802e-13,-11450.6,2.70709], Tmin=(782.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.3373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsJOOC)"""),
)

species(
    label = 'COOC=[C]OOF(3405)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {1,S} {4,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {8,D} {12,S}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (208.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.896451,'amu*angstrom^2'), symmetry=1, barrier=(20.6112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894899,'amu*angstrom^2'), symmetry=1, barrier=(20.5755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.737738,0.0778002,-0.000111671,8.87863e-08,-2.87527e-11,25191.3,31.6473], Tmin=(100,'K'), Tmax=(752.117,'K')), NASAPolynomial(coeffs=[9.86634,0.0292435,-1.48148e-05,2.92022e-09,-2.06521e-13,23818.4,-9.79099], Tmin=(752.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
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
    E0 = (-109.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (277.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-56.5548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-132.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-61.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-15.3413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-103.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (329.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-134.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (146.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (148.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (34.4668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (263.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (309.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (312.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-9.52879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-58.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (292.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['COOC=C(F)O[O](3351)'],
    products = ['CH2O(20)', '[O]OC(F)C=O(1537)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1054.45,'s^-1'), n=2.72887, Ea=(31.87,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[O]OC(F)=COO(3391)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.67448e-06,'m^3/(mol*s)'), n=3.15288, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node OH_2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['COOC=C(F)O[O](3351)'],
    products = ['COC(F)(C=O)O[O](3392)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(85.022,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'COO[CH]C(=O)F(3393)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['COOC=C(F)O[O](3351)'],
    products = ['CH3O(27)', 'FC1=COOO1(3394)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(80.3459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 78.4 to 80.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['COOC=C(F)O[O](3351)'],
    products = ['COOC1OO[C]1F(3395)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['COOC=C(F)O[O](3351)'],
    products = ['COO[CH]C1(F)OO1(3396)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(38.4928,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'COOC=[C]O[O](3397)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3O(27)', '[O]O[C](F)C=O(1584)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(52.8868,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(19)', '[O]OC=C(F)O[O](3398)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(284464,'m^3/(mol*s)'), n=0.158989, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.07307342555193728, var=1.039002618272806, Tref=1000.0, N=14, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO[O](183)', '[CH]=C(F)O[O](1081)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O2(2)', 'COOC=[C]F(3399)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[CH2]OOC=C(F)O[O](3400)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'COO[C]=C(F)O[O](3401)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'COOC#CO[O](3402)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(286.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction16',
    reactants = ['COO[C]=C(F)OO(3403)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['COOC=C(F)O[O](3351)'],
    products = ['[CH2]OOC=C(F)OO(3404)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(271800,'s^-1'), n=1.51, Ea=(83.4708,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 266 used for R7H_OOCs4;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R7H_OOCs4;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['COOC=[C]OOF(3405)'],
    products = ['COOC=C(F)O[O](3351)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.8166,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1174',
    isomers = [
        'COOC=C(F)O[O](3351)',
    ],
    reactants = [
        ('CH2O(20)', '[O]OC(F)C=O(1537)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1174',
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

