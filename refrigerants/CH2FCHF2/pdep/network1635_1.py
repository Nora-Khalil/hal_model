species(
    label = '[O]OC(F)C(F)OOC(F)=CF(5351)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {6,S} {12,D}
12 C u0 p0 c0 {4,S} {11,D} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-804.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.253806,'amu*angstrom^2'), symmetry=1, barrier=(5.8355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83275,'amu*angstrom^2'), symmetry=1, barrier=(42.1386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83382,'amu*angstrom^2'), symmetry=1, barrier=(42.1632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252969,'amu*angstrom^2'), symmetry=1, barrier=(5.81625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83436,'amu*angstrom^2'), symmetry=1, barrier=(42.1755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22385,0.126677,-0.000208037,1.8417e-07,-6.51911e-11,-96631.8,39.2086], Tmin=(100,'K'), Tmax=(765.687,'K')), NASAPolynomial(coeffs=[12.3877,0.0442385,-2.43391e-05,4.90009e-09,-3.48259e-13,-98384,-20.6562], Tmin=(765.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-804.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C(=O)F(4704)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-582.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.519617,'amu*angstrom^2'), symmetry=1, barrier=(11.947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54035,'amu*angstrom^2'), symmetry=1, barrier=(35.4156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3788.12,'J/mol'), sigma=(5.7366,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=591.70 K, Pc=45.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66217,0.0558669,-8.74757e-05,7.36244e-08,-2.48646e-11,-69953.7,20.9883], Tmin=(100,'K'), Tmax=(757.994,'K')), NASAPolynomial(coeffs=[8.40288,0.0184184,-9.65371e-06,1.91163e-09,-1.34843e-13,-70921.6,-9.30828], Tmin=(757.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ)"""),
)

species(
    label = 'O=C(F)CF(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-611.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,1102,1254,1379,1468,1475,3037,3083,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.393549,'amu*angstrom^2'), symmetry=1, barrier=(15.7349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84069,0.0198378,-9.07231e-06,-3.50907e-10,9.31991e-13,-73601.2,9.53853], Tmin=(10,'K'), Tmax=(1222.46,'K')), NASAPolynomial(coeffs=[7.66923,0.0123956,-6.18005e-06,1.47455e-09,-1.37207e-13,-74917.2,-11.2551], Tmin=(1222.46,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-611.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OC(F)OOC(F)=CF(6079)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {3,S} {9,D} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-598.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,509,613,660,1171,1360,1414,3084,326,540,652,719,1357,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.409887,'amu*angstrom^2'), symmetry=1, barrier=(9.42411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410886,'amu*angstrom^2'), symmetry=1, barrier=(9.44708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.304827,0.107689,-0.000194569,1.8062e-07,-6.44937e-11,-71850.5,33.0944], Tmin=(100,'K'), Tmax=(827.449,'K')), NASAPolynomial(coeffs=[9.48745,0.0358036,-1.97542e-05,3.91857e-09,-2.73651e-13,-72630.6,-7.21486], Tmin=(827.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = 'HO2(11)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.49012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1064.4,1465.7,3224.93],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'F[C]C(F)OOC(F)=CF(3339)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {12,S}
10 C u0 p1 c0 {4,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-479.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,326,540,652,719,1357,194,682,905,1196,1383,3221,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.59354,'amu*angstrom^2'), symmetry=1, barrier=(36.6386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.515995,'amu*angstrom^2'), symmetry=1, barrier=(11.8637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59264,'amu*angstrom^2'), symmetry=1, barrier=(36.618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59349,'amu*angstrom^2'), symmetry=1, barrier=(36.6375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0399941,0.0990098,-0.000167034,1.51856e-07,-5.47798e-11,-57583.6,32.4423], Tmin=(100,'K'), Tmax=(779.846,'K')), NASAPolynomial(coeffs=[9.59859,0.0360867,-2.00669e-05,4.04579e-09,-2.8731e-13,-58676.9,-9.0326], Tmin=(779.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]OC(OOC(F)=CF)C(F)F(6080)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {9,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {5,S} {7,S} {10,S} {13,S}
10 C u0 p0 c0 {1,S} {2,S} {9,S} {14,S}
11 C u0 p0 c0 {3,S} {6,S} {12,D}
12 C u0 p0 c0 {4,S} {11,D} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-818.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3828,0.133806,-0.00023718,2.1896e-07,-7.83553e-11,-98313.3,41.2479], Tmin=(100,'K'), Tmax=(816.805,'K')), NASAPolynomial(coeffs=[11.0446,0.0456544,-2.51735e-05,5.01369e-09,-3.5172e-13,-99433,-10.619], Tmin=(816.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-818.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C(F)OC(F)C(=O)F(6081)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {12,D}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {5,S} {12,S} {15,S}
12 C u0 p0 c0 {4,S} {7,D} {11,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1195.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.802264,0.115283,-0.000161901,1.19744e-07,-3.57978e-11,-143641,36.4413], Tmin=(100,'K'), Tmax=(813.447,'K')), NASAPolynomial(coeffs=[14.736,0.0388751,-2.10033e-05,4.2682e-09,-3.07923e-13,-146169,-35.3135], Tmin=(813.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1195.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ)"""),
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
    label = '[O]C(F)C(F)OOC(F)=CF(6082)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {4,S} {10,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-801.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,391,562,707,872,1109,1210,1289,3137,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.190825,'amu*angstrom^2'), symmetry=1, barrier=(4.38744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14622,'amu*angstrom^2'), symmetry=1, barrier=(49.3458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190711,'amu*angstrom^2'), symmetry=1, barrier=(4.38482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.95701,'amu*angstrom^2'), symmetry=1, barrier=(113.971,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (175.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.505073,0.108636,-0.00016982,1.46461e-07,-5.13228e-11,-96286.9,36.0151], Tmin=(100,'K'), Tmax=(731.925,'K')), NASAPolynomial(coeffs=[11.2811,0.0403265,-2.18381e-05,4.39773e-09,-3.13707e-13,-97907.8,-16.455], Tmin=(731.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-801.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'O=C(F)[CH]F(215)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-442.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(1.5365,'amu*angstrom^2'), symmetry=1, barrier=(35.3272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3038.52,'J/mol'), sigma=(4.81134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=474.61 K, Pc=61.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91969,0.0052474,7.61145e-05,-1.81228e-07,1.26791e-10,-53179.7,10.2212], Tmin=(10,'K'), Tmax=(489.215,'K')), NASAPolynomial(coeffs=[4.0714,0.0192239,-1.3397e-05,4.3331e-09,-5.27008e-13,-53376.6,7.73666], Tmin=(489.215,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-442.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1OOOC1F(909)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-475.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91481,0.0309223,2.38538e-05,-6.78047e-08,3.32829e-11,-57125.8,14.9275], Tmin=(100,'K'), Tmax=(890.994,'K')), NASAPolynomial(coeffs=[17.3157,0.00257128,2.91375e-06,-7.56452e-10,5.35431e-14,-61489.3,-66.6814], Tmin=(890.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-475.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(123trioxolane)"""),
)

species(
    label = 'FC=C(F)OOC(F)=CF(4776)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u0 p0 c0 {2,S} {6,S} {10,D}
9  C u0 p0 c0 {3,S} {7,D} {11,S}
10 C u0 p0 c0 {4,S} {8,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-641.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,262,390,483,597,572,732,631,807,1275,1439,151,237,609,755,844,966,1147,1245,1323,1443,3181,3261,180],'cm^-1')),
        HinderedRotor(inertia=(0.138574,'amu*angstrom^2'), symmetry=1, barrier=(3.18609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138096,'amu*angstrom^2'), symmetry=1, barrier=(3.17509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13868,'amu*angstrom^2'), symmetry=1, barrier=(3.18854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.143224,0.0993805,-0.000183518,1.7913e-07,-6.68155e-11,-76970.6,32.2265], Tmin=(100,'K'), Tmax=(824.796,'K')), NASAPolynomial(coeffs=[5.75632,0.0416549,-2.3061e-05,4.59602e-09,-3.22245e-13,-76859,12.5177], Tmin=(824.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-641.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'F[C]1OOC(F)C(F)OOC1F(6083)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {6,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {7,S} {12,S} {15,S}
12 C u1 p0 c0 {4,S} {8,S} {11,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-863.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.413761,0.0825041,-5.44324e-05,1.13469e-08,5.33085e-13,-103680,30.4199], Tmin=(100,'K'), Tmax=(1285.06,'K')), NASAPolynomial(coeffs=[23.3819,0.0285493,-1.4931e-05,3.0342e-09,-2.19246e-13,-111457,-96.8102], Tmin=(1285.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-863.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(CsCFHO) + group(CsCFHO) + group(CsCFHO) + ring(Cyclooctane) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C1(F)OOC(F)C(F)OO1(6084)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {7,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
12 C u1 p0 c0 {4,S} {11,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-902.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30213,0.0729591,5.01081e-05,-1.60684e-07,7.82068e-11,-108267,33.4745], Tmin=(100,'K'), Tmax=(927.069,'K')), NASAPolynomial(coeffs=[45.4494,-0.0161214,1.19933e-05,-2.21866e-09,1.3233e-13,-121776,-214.64], Tmin=(927.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-902.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cycloheptane) + radical(Csj(Cs-F1sO2sO2s)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = '[O]O[CH]C(F)OOC(F)=CF(6085)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {6,S} {8,S} {13,S}
10 C u0 p0 c0 {2,S} {5,S} {11,D}
11 C u0 p0 c0 {3,S} {10,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-405.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,326,540,652,719,1357,194,682,905,1196,1383,3221,180,2634.51],'cm^-1')),
        HinderedRotor(inertia=(2.0959,'amu*angstrom^2'), symmetry=1, barrier=(48.1888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401462,'amu*angstrom^2'), symmetry=1, barrier=(9.23039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09971,'amu*angstrom^2'), symmetry=1, barrier=(48.2764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400839,'amu*angstrom^2'), symmetry=1, barrier=(9.21607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09026,'amu*angstrom^2'), symmetry=1, barrier=(48.0591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.96551,0.125458,-0.000228203,2.14982e-07,-7.76287e-11,-48582.4,37.8888], Tmin=(100,'K'), Tmax=(831.079,'K')), NASAPolynomial(coeffs=[9.02564,0.0450276,-2.46592e-05,4.87774e-09,-3.40053e-13,-49126.1,-1.74414], Tmin=(831.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-405.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[O]OC(F)[CH]OOC(F)=CF(6086)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
9  C u1 p0 c0 {4,S} {8,S} {13,S}
10 C u0 p0 c0 {2,S} {5,S} {11,D}
11 C u0 p0 c0 {3,S} {10,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-407.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.235541,'amu*angstrom^2'), symmetry=1, barrier=(5.41555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30824,'amu*angstrom^2'), symmetry=1, barrier=(53.071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237753,'amu*angstrom^2'), symmetry=1, barrier=(5.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235761,'amu*angstrom^2'), symmetry=1, barrier=(5.4206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.3036,'amu*angstrom^2'), symmetry=1, barrier=(52.9643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.929875,0.124189,-0.000224184,2.10498e-07,-7.5949e-11,-48840.8,37.9068], Tmin=(100,'K'), Tmax=(828.73,'K')), NASAPolynomial(coeffs=[9.13033,0.0447353,-2.44512e-05,4.83813e-09,-3.37566e-13,-49447.3,-2.33681], Tmin=(828.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CCsJOOC)"""),
)

species(
    label = '[O]OC(F)C(F)OO[C]=CF(6087)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {11,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {11,D} {14,S}
11 C u1 p0 c0 {6,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-378.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,615,860,1140,1343,3152,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702294,0.113288,-0.000171242,1.38281e-07,-4.51614e-11,-45406.1,36.0518], Tmin=(100,'K'), Tmax=(746.749,'K')), NASAPolynomial(coeffs=[13.2512,0.0385533,-2.1139e-05,4.28915e-09,-3.07926e-13,-47490.3,-27.1925], Tmin=(746.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-378.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(Cdj(Cd-F1sH)(O2s-O2s))"""),
)

species(
    label = '[CH]=C(F)OOC(F)C(F)O[O](6088)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {5,S} {11,D}
11 C u1 p0 c0 {10,D} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-375.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,293,496,537,1218,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.186045,'amu*angstrom^2'), symmetry=1, barrier=(4.27755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.573832,'amu*angstrom^2'), symmetry=1, barrier=(13.1935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81367,'amu*angstrom^2'), symmetry=1, barrier=(41.6998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81493,'amu*angstrom^2'), symmetry=1, barrier=(41.7288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81554,'amu*angstrom^2'), symmetry=1, barrier=(41.7428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.901642,0.118109,-0.000191089,1.65298e-07,-5.74068e-11,-45006.4,37.3964], Tmin=(100,'K'), Tmax=(745.461,'K')), NASAPolynomial(coeffs=[12.8787,0.0389986,-2.15072e-05,4.34119e-09,-3.09347e-13,-46917.3,-24.0746], Tmin=(745.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-375.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = '[O]OC(F)C([O])F(4580)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u1 p2 c0 {7,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-403.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(0.434449,'amu*angstrom^2'), symmetry=1, barrier=(9.98883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.938918,'amu*angstrom^2'), symmetry=1, barrier=(21.5876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61323,0.0562001,-7.44595e-05,5.23724e-08,-1.48508e-11,-48394.3,22.2333], Tmin=(100,'K'), Tmax=(857.912,'K')), NASAPolynomial(coeffs=[9.58518,0.0190297,-9.46747e-06,1.86664e-09,-1.3269e-13,-49762.1,-15.0048], Tmin=(857.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(O2sj(Cs-CsF1sH)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=CF(214)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {2,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-249.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.358734,'amu*angstrom^2'), symmetry=1, barrier=(8.24799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3441.78,'J/mol'), sigma=(5.35945,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.60 K, Pc=50.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75962,0.0203814,3.21936e-05,-1.12324e-07,8.28e-11,-30007.8,12.0802], Tmin=(10,'K'), Tmax=(521.492,'K')), NASAPolynomial(coeffs=[5.63612,0.021414,-1.5147e-05,4.91851e-09,-5.9756e-13,-30413.3,2.2378], Tmin=(521.492,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-249.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)[CH]F(131)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-229.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.418102,'amu*angstrom^2'), symmetry=1, barrier=(9.61299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419586,'amu*angstrom^2'), symmetry=1, barrier=(9.64711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.65,'J/mol'), sigma=(5.51737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.71 K, Pc=45.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66565,0.0577513,-0.000101954,9.14652e-08,-3.12262e-11,-27473.1,20.2391], Tmin=(100,'K'), Tmax=(863.112,'K')), NASAPolynomial(coeffs=[7.65508,0.0169238,-8.28535e-06,1.57109e-09,-1.06542e-13,-28020.1,-4.95442], Tmin=(863.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'CHFCF[Z](72)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-58.4765,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([18.6764,93.2568,111.933],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([203.406,408.59,750.559,801.84,1025.92,1150.34,1361.98,1782.01,3238.53],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9644,0.00217073,4.76141e-05,-9.63411e-08,5.92621e-11,-7031.25,9.1333], Tmin=(10,'K'), Tmax=(529.917,'K')), NASAPolynomial(coeffs=[3.26037,0.0150322,-1.01557e-05,3.21347e-09,-3.84694e-13,-7062.61,11.0829], Tmin=(529.917,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-58.4765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""F[C]DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)C(F)O[O](6089)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {3,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-406.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,180],'cm^-1')),
        HinderedRotor(inertia=(0.732029,'amu*angstrom^2'), symmetry=1, barrier=(16.8308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.73179,'amu*angstrom^2'), symmetry=1, barrier=(16.8253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13009,'amu*angstrom^2'), symmetry=1, barrier=(48.9749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866409,0.0745832,-0.000113958,9.19689e-08,-2.96949e-11,-48738,24.8336], Tmin=(100,'K'), Tmax=(758.48,'K')), NASAPolynomial(coeffs=[10.64,0.0230422,-1.20321e-05,2.38504e-09,-1.68639e-13,-50220.7,-19.617], Tmin=(758.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(ROOJ) + radical(ROOJ)"""),
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
    label = 'F[CH]C(F)OOC(F)=CF(3348)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u1 p0 c0 {3,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-633.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06785,'amu*angstrom^2'), symmetry=1, barrier=(47.5439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452614,'amu*angstrom^2'), symmetry=1, barrier=(10.4065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06733,'amu*angstrom^2'), symmetry=1, barrier=(47.532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245044,'amu*angstrom^2'), symmetry=1, barrier=(5.63404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.059,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.85,'J/mol'), sigma=(5.29792,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=518.55 K, Pc=50.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.382469,0.106753,-0.000177615,1.58956e-07,-5.65336e-11,-75994.1,33.0843], Tmin=(100,'K'), Tmax=(780.036,'K')), NASAPolynomial(coeffs=[10.6309,0.0377033,-2.06549e-05,4.14348e-09,-2.93527e-13,-77329.7,-14.8611], Tmin=(780.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-633.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[O]O[CH]F(209)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-10.1775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,580,1155,1237,1373,3147],'cm^-1')),
        HinderedRotor(inertia=(0.627597,'amu*angstrom^2'), symmetry=1, barrier=(14.4297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87623,0.0275715,-4.75678e-05,4.25822e-08,-1.46837e-11,-1186.36,12.8042], Tmin=(100,'K'), Tmax=(836.137,'K')), NASAPolynomial(coeffs=[5.84968,0.00835189,-4.1278e-06,8.02249e-10,-5.55883e-14,-1508.99,0.0352551], Tmin=(836.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.1775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'F[CH]OOC(F)=CF(3245)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {9,S}
8  C u1 p0 c0 {3,S} {5,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-409.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,326,540,652,719,1357,194,682,905,1196,1383,3221,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.226767,'amu*angstrom^2'), symmetry=1, barrier=(5.21382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22726,'amu*angstrom^2'), symmetry=1, barrier=(5.22516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.46124,'amu*angstrom^2'), symmetry=1, barrier=(56.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91761,0.0785622,-0.000142421,1.3701e-07,-5.06695e-11,-49160.5,26.5988], Tmin=(100,'K'), Tmax=(821.716,'K')), NASAPolynomial(coeffs=[6.02591,0.0319352,-1.75834e-05,3.5008e-09,-2.45504e-13,-49265.4,7.42757], Tmin=(821.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-409.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFHHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsF1sHO2s)"""),
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
    label = '[O]OC(F)[C](F)OOC(F)=CF(6090)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {9,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {5,S} {9,S}
11 C u0 p0 c0 {3,S} {6,S} {12,D}
12 C u0 p0 c0 {4,S} {11,D} {14,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-610.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,395,473,707,1436,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.20835,'amu*angstrom^2'), symmetry=1, barrier=(50.7742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376941,'amu*angstrom^2'), symmetry=1, barrier=(8.66662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.377589,'amu*angstrom^2'), symmetry=1, barrier=(8.68151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37774,'amu*angstrom^2'), symmetry=1, barrier=(8.68498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21347,'amu*angstrom^2'), symmetry=1, barrier=(50.8921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38977,0.135595,-0.000249767,2.33787e-07,-8.38218e-11,-73231.5,40.6044], Tmin=(100,'K'), Tmax=(828.814,'K')), NASAPolynomial(coeffs=[10.7282,0.0442062,-2.48172e-05,4.944e-09,-3.45674e-13,-74110,-8.76458], Tmin=(828.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]O[C](F)C(F)OOC(F)=CF(6091)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {7,S} {9,S}
11 C u0 p0 c0 {3,S} {6,S} {12,D}
12 C u0 p0 c0 {4,S} {11,D} {14,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-610.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,395,473,707,1436,326,540,652,719,1357,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.20835,'amu*angstrom^2'), symmetry=1, barrier=(50.7742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376941,'amu*angstrom^2'), symmetry=1, barrier=(8.66662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.377589,'amu*angstrom^2'), symmetry=1, barrier=(8.68151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37774,'amu*angstrom^2'), symmetry=1, barrier=(8.68498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21347,'amu*angstrom^2'), symmetry=1, barrier=(50.8921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38977,0.135595,-0.000249767,2.33787e-07,-8.38218e-11,-73231.5,40.6044], Tmin=(100,'K'), Tmax=(828.814,'K')), NASAPolynomial(coeffs=[10.7282,0.0442062,-2.48172e-05,4.944e-09,-3.45674e-13,-74110,-8.76458], Tmin=(828.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)C(F)OOC(F)=[C]F(6092)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {6,S} {12,D}
12 C u1 p0 c0 {4,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-535.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,293,496,537,1218,167,640,1190,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22918,0.129043,-0.000226134,2.07399e-07,-7.44069e-11,-64214.8,40.2618], Tmin=(100,'K'), Tmax=(799.071,'K')), NASAPolynomial(coeffs=[11.4923,0.0434961,-2.45019e-05,4.93256e-09,-3.4879e-13,-65549.8,-13.8905], Tmin=(799.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-535.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
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
    label = '[O]OC=COOC(F)=CF(6093)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {3,S} {10,D}
9  C u0 p0 c0 {5,S} {7,D} {12,S}
10 C u0 p0 c0 {2,S} {8,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-221.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.438704,'amu*angstrom^2'), symmetry=1, barrier=(10.0867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438298,'amu*angstrom^2'), symmetry=1, barrier=(10.0773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.73721,'amu*angstrom^2'), symmetry=1, barrier=(16.9499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438334,'amu*angstrom^2'), symmetry=1, barrier=(10.0782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.14715,0.09817,-0.00014957,1.19823e-07,-3.82415e-11,-26464.3,35.7625], Tmin=(100,'K'), Tmax=(768.601,'K')), NASAPolynomial(coeffs=[13.1193,0.0291172,-1.47861e-05,2.89701e-09,-2.03396e-13,-28503.3,-24.747], Tmin=(768.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
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
    label = '[O]OC=C(F)OOC(F)=CF(6094)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {2,S} {5,S} {10,D}
9  C u0 p0 c0 {1,S} {4,S} {11,D}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 C u0 p0 c0 {3,S} {9,D} {13,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-379.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.219995,'amu*angstrom^2'), symmetry=1, barrier=(5.05812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220404,'amu*angstrom^2'), symmetry=1, barrier=(5.06753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219927,'amu*angstrom^2'), symmetry=1, barrier=(5.05656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220388,'amu*angstrom^2'), symmetry=1, barrier=(5.06715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.375486,0.114773,-0.000221195,2.18226e-07,-8.10014e-11,-45520,38.585], Tmin=(100,'K'), Tmax=(840.733,'K')), NASAPolynomial(coeffs=[5.4139,0.046189,-2.56085e-05,5.07154e-09,-3.52898e-13,-45043.1,20.2848], Tmin=(840.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=COOC(F)=CF(4950)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {5,S} {10,D} {12,S}
9  C u0 p0 c0 {1,S} {4,S} {11,D}
10 C u0 p0 c0 {2,S} {6,S} {8,D}
11 C u0 p0 c0 {3,S} {9,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-379.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,262,390,483,597,572,732,631,807,1275,1439,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.219995,'amu*angstrom^2'), symmetry=1, barrier=(5.05812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220404,'amu*angstrom^2'), symmetry=1, barrier=(5.06753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219927,'amu*angstrom^2'), symmetry=1, barrier=(5.05656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220388,'amu*angstrom^2'), symmetry=1, barrier=(5.06715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.375486,0.114773,-0.000221195,2.18226e-07,-8.10014e-11,-45520,38.585], Tmin=(100,'K'), Tmax=(840.733,'K')), NASAPolynomial(coeffs=[5.4139,0.046189,-2.56085e-05,5.07154e-09,-3.52898e-13,-45043.1,20.2848], Tmin=(840.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = 'C#COOC(F)C(F)O[O](6095)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {3,S} {9,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {10,T}
10 C u0 p0 c0 {9,T} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-221.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.53251,'amu*angstrom^2'), symmetry=1, barrier=(35.2355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53265,'amu*angstrom^2'), symmetry=1, barrier=(35.2385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53206,'amu*angstrom^2'), symmetry=1, barrier=(35.2251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53203,'amu*angstrom^2'), symmetry=1, barrier=(35.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53212,'amu*angstrom^2'), symmetry=1, barrier=(35.2265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.265205,0.101776,-0.000141348,1.01159e-07,-2.90508e-11,-26471.9,31.2239], Tmin=(100,'K'), Tmax=(847.855,'K')), NASAPolynomial(coeffs=[14.617,0.031565,-1.7134e-05,3.48998e-09,-2.52298e-13,-28995.5,-38.1182], Tmin=(847.855,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C(F)OOC#CF(6096)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
9  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {6,S} {11,T}
11 C u0 p0 c0 {3,S} {10,T}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-321.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,2175,525,239,401,1367,180,304.328],'cm^-1')),
        HinderedRotor(inertia=(1.59562,'amu*angstrom^2'), symmetry=1, barrier=(36.6865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60111,'amu*angstrom^2'), symmetry=1, barrier=(36.8126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651435,'amu*angstrom^2'), symmetry=1, barrier=(14.9778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59935,'amu*angstrom^2'), symmetry=1, barrier=(36.7721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59883,'amu*angstrom^2'), symmetry=1, barrier=(36.7603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.714421,0.113754,-0.000182737,1.54913e-07,-5.28016e-11,-38457.3,35.1396], Tmin=(100,'K'), Tmax=(717.731,'K')), NASAPolynomial(coeffs=[13.2215,0.0360923,-2.04413e-05,4.1743e-09,-2.99756e-13,-40457.9,-27.4723], Tmin=(717.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtOs) + group(CtCF) + radical(ROOJ)"""),
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
    label = 'O=C(F)C(F)OOC(F)=CF(6097)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {9,D}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {4,S} {10,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-981.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395264,0.107514,-0.000179662,1.62782e-07,-5.87337e-11,-117849,34.556], Tmin=(100,'K'), Tmax=(774.637,'K')), NASAPolynomial(coeffs=[10.1024,0.0397043,-2.20159e-05,4.4403e-09,-3.15631e-13,-119067,-10.7741], Tmin=(774.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-981.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'OOC(F)[C](F)OOC(F)=CF(6098)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {6,S} {11,S}
8  O u0 p2 c0 {5,S} {15,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {6,S} {9,S}
11 C u0 p0 c0 {3,S} {7,S} {12,D}
12 C u0 p0 c0 {4,S} {11,D} {14,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-762.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,500,795,815,487,638,688,1119,1325,1387,3149,395,473,707,1436,326,540,652,719,1357,194,682,905,1196,1383,3221,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62707,0.138888,-0.000244168,2.23744e-07,-8.00606e-11,-91503,40.6968], Tmin=(100,'K'), Tmax=(801.99,'K')), NASAPolynomial(coeffs=[12.1926,0.0460855,-2.59378e-05,5.21592e-09,-3.68426e-13,-92951.8,-18.1385], Tmin=(801.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-762.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OOC(F)C(F)OOC(F)=[C]F(6099)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {6,S} {15,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {7,S} {12,D}
12 C u1 p0 c0 {4,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-687.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,293,496,537,1218,167,640,1190,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42599,0.131823,-0.000218602,1.94733e-07,-6.95356e-11,-82488,40.2115], Tmin=(100,'K'), Tmax=(752.212,'K')), NASAPolynomial(coeffs=[12.7963,0.0456655,-2.57973e-05,5.24713e-09,-3.75165e-13,-84329.8,-22.3735], Tmin=(752.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'FC=C(F)OO[CH]C(F)OOF(6100)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {6,S} {11,S}
8  O u0 p2 c0 {4,S} {5,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {6,S} {9,S} {14,S}
11 C u0 p0 c0 {2,S} {7,S} {12,D}
12 C u0 p0 c0 {3,S} {11,D} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-464.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,500,795,815,277,555,632,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(1.74066,'amu*angstrom^2'), symmetry=1, barrier=(40.0213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73968,'amu*angstrom^2'), symmetry=1, barrier=(39.9987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74307,'amu*angstrom^2'), symmetry=1, barrier=(40.0765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664123,'amu*angstrom^2'), symmetry=1, barrier=(15.2695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73815,'amu*angstrom^2'), symmetry=1, barrier=(39.9636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53527,0.137257,-0.000241764,2.2332e-07,-8.0422e-11,-55638.7,40.7012], Tmin=(100,'K'), Tmax=(805.641,'K')), NASAPolynomial(coeffs=[11.2807,0.0476774,-2.66661e-05,5.34883e-09,-3.77239e-13,-56861.6,-13.1324], Tmin=(805.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-464.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCsJOOC)"""),
)

species(
    label = 'FC=C(F)OOC(F)[CH]OOF(6101)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {4,S} {7,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u0 p0 c0 {2,S} {6,S} {12,D}
12 C u0 p0 c0 {3,S} {11,D} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-463.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.674982,'amu*angstrom^2'), symmetry=1, barrier=(15.5192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72979,'amu*angstrom^2'), symmetry=1, barrier=(39.7714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72748,'amu*angstrom^2'), symmetry=1, barrier=(39.7181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72844,'amu*angstrom^2'), symmetry=1, barrier=(39.7403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73162,'amu*angstrom^2'), symmetry=1, barrier=(39.8134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55445,0.13792,-0.000243888,2.25739e-07,-8.13472e-11,-55517.2,40.6976], Tmin=(100,'K'), Tmax=(807.19,'K')), NASAPolynomial(coeffs=[11.2231,0.0478472,-2.67894e-05,5.37283e-09,-3.78783e-13,-56708.4,-12.8111], Tmin=(807.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-463.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCsJOO)"""),
)

species(
    label = 'FC=[C]OOC(F)C(F)OOF(6102)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {12,S}
8  O u0 p2 c0 {4,S} {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {12,D} {15,S}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-435.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,615,860,1140,1343,3152,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.59763,'amu*angstrom^2'), symmetry=1, barrier=(36.7326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5978,'amu*angstrom^2'), symmetry=1, barrier=(36.7365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59717,'amu*angstrom^2'), symmetry=1, barrier=(36.7221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59621,'amu*angstrom^2'), symmetry=1, barrier=(36.7,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59752,'amu*angstrom^2'), symmetry=1, barrier=(36.7301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59757,'amu*angstrom^2'), symmetry=1, barrier=(36.7314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1062,0.123702,-0.000178066,1.34629e-07,-4.12501e-11,-52212.6,38.1404], Tmin=(100,'K'), Tmax=(793.156,'K')), NASAPolynomial(coeffs=[15.0775,0.0420838,-2.37091e-05,4.88655e-09,-3.54967e-13,-54779.8,-36.1859], Tmin=(793.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-F1sH)(O2s-O2s))"""),
)

species(
    label = '[CH]=C(F)OOC(F)C(F)OOF(6103)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {4,S} {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {7,S} {12,D}
12 C u1 p0 c0 {11,D} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-432.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,293,496,537,1218,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.51792,'amu*angstrom^2'), symmetry=1, barrier=(34.9,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51768,'amu*angstrom^2'), symmetry=1, barrier=(34.8945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51801,'amu*angstrom^2'), symmetry=1, barrier=(34.9021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51814,'amu*angstrom^2'), symmetry=1, barrier=(34.9049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51771,'amu*angstrom^2'), symmetry=1, barrier=(34.8951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30789,0.12849,-0.000197467,1.60513e-07,-5.27181e-11,-51812.7,39.4968], Tmin=(100,'K'), Tmax=(742.735,'K')), NASAPolynomial(coeffs=[14.6228,0.0426903,-2.418e-05,4.96462e-09,-3.58662e-13,-54179.1,-32.6206], Tmin=(742.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-432.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sO2s)(H))"""),
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
    E0 = (-437.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (22.4547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-80.2365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-156.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-373.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-209.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-400.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-284.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-397.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-359.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (16.8868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (14.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (43.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (46.5895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-446.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-129.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-115.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-292.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-70.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-49.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-49.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (25.7282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (119.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-114.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-114.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (119.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (18.7035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-298.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-329.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-189.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-45.4966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-21.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-23.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-21.5192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['[O]OC(F)C(=O)F(4704)', 'O=C(F)CF(879)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(18.122,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[O]OC(F)OOC(F)=CF(6079)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.25106e-07,'m^3/(mol*s)'), n=3.34134, Ea=(132.994,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['HO2(11)', 'F[C]C(F)OOC(F)=CF(3339)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(47.8902,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['[O]OC(OOC(F)=CF)C(F)F(6080)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['[O]OC(F)C(F)OC(F)C(=O)F(6081)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(82.5351,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', '[O]C(F)C(F)OOC(F)=CF(6082)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['O=C(F)[CH]F(215)', 'FC1OOOC1F(909)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['HO2(11)', 'FC=C(F)OOC(F)=CF(4776)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(170.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['F[C]1OOC(F)C(F)OOC1F(6083)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.97015e+10,'s^-1'), n=0.301339, Ea=(58.2388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra;radadd_intra] for rate rule [R8_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['F[CH]C1(F)OOC(F)C(F)OO1(6084)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(95.8805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;doublebond_intra;radadd_intra_O] + [R8;doublebond_intra;radadd_intra] for rate rule [R8;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]O[CH]C(F)OOC(F)=CF(6085)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]OC(F)[CH]OOC(F)=CF(6086)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]OC(F)C(F)OO[C]=CF(6087)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[CH]=C(F)OOC(F)C(F)O[O](6088)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C(F)[CH]F(215)', '[O]OC(F)C([O])F(4580)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(49.7802,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC(F)=CF(214)', '[O]OC(F)[CH]F(131)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHFCF[Z](72)', '[O]OC(F)C(F)O[O](6089)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O2(2)', 'F[CH]C(F)OOC(F)=CF(3348)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]O[CH]F(209)', 'F[CH]OOC(F)=CF(3245)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]OC(F)[C](F)OOC(F)=CF(6090)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]O[C](F)C(F)OOC(F)=CF(6091)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', '[O]OC(F)C(F)OOC(F)=[C]F(6092)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F2(78)', '[O]OC=COOC(F)=CF(6093)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', '[O]OC=C(F)OOC(F)=CF(6094)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(196.849,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', '[O]OC(F)=COOC(F)=CF(4950)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(196.849,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F2(78)', 'C#COOC(F)C(F)O[O](6095)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', '[O]OC(F)C(F)OOC#CF(6096)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(271.613,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    products = ['OH(7)', 'O=C(F)C(F)OOC(F)=CF(6097)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['OOC(F)[C](F)OOC(F)=CF(6098)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;O_H_out] for rate rule [R4H_SSS;C_rad_out_noH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['OOC(F)C(F)OOC(F)=[C]F(6099)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.20037e+07,'s^-1'), n=1.48577, Ea=(148.291,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_single;XH_out] + [R8H;Y_rad_out;XH_out] for rate rule [R8H;Cd_rad_out_single;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['FC=C(F)OO[CH]C(F)OOF(6100)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(69.3539,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction32',
    reactants = ['FC=C(F)OOC(F)[CH]OOF(6101)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.4141,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FC=[C]OOC(F)C(F)OOF(6102)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(62.3436,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(F)OOC(F)C(F)OOF(6103)'],
    products = ['[O]OC(F)C(F)OOC(F)=CF(5351)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(61.4869,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1635',
    isomers = [
        '[O]OC(F)C(F)OOC(F)=CF(5351)',
    ],
    reactants = [
        ('[O]OC(F)C(=O)F(4704)', 'O=C(F)CF(879)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1635',
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

