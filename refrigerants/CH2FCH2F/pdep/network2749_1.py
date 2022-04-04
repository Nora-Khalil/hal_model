species(
    label = 'FC=[C]C(F)(F)C[C](F)F(7492)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-610.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,136,307,446,511,682,757,1180,1185,190,488,555,1236,1407,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.214946,'amu*angstrom^2'), symmetry=1, barrier=(4.94204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21419,'amu*angstrom^2'), symmetry=1, barrier=(4.92464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212168,'amu*angstrom^2'), symmetry=1, barrier=(4.87815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.192983,0.103232,-0.000175923,1.60091e-07,-5.70865e-11,-73254.7,32.6261], Tmin=(100,'K'), Tmax=(807.521,'K')), NASAPolynomial(coeffs=[9.55458,0.037105,-1.99452e-05,3.95658e-09,-2.77761e-13,-74247.2,-8.71416], Tmin=(807.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'CH2CF2(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.616,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=C=C(F)F(5948)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-364.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,94,120,354,641,825,1294,540,610,2055,2775.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.43,'J/mol'), sigma=(4.61545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.73 K, Pc=67.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78631,0.0191407,1.84311e-05,-5.99612e-08,3.72113e-11,-43861.2,10.8341], Tmin=(10,'K'), Tmax=(620.754,'K')), NASAPolynomial(coeffs=[5.19832,0.0213405,-1.41861e-05,4.38954e-09,-5.1374e-13,-44254.2,2.94178], Tmin=(620.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-364.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FCDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)(F)C(F)(F)[C]=CF(7489)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
8  C u1 p0 c0 {6,S} {11,S} {12,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-600.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,136,307,446,511,682,757,1180,1185,3000,3100,440,815,1455,1000,615,860,1140,1343,3152,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.176693,'amu*angstrom^2'), symmetry=1, barrier=(4.06252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176407,'amu*angstrom^2'), symmetry=1, barrier=(4.05595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614109,'amu*angstrom^2'), symmetry=1, barrier=(14.1196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377058,0.104991,-0.000170601,1.46261e-07,-4.97548e-11,-72104.2,32.3579], Tmin=(100,'K'), Tmax=(782.485,'K')), NASAPolynomial(coeffs=[12.3823,0.0323484,-1.71297e-05,3.39234e-09,-2.38619e-13,-73873.9,-24.6183], Tmin=(782.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-600.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'F[CH]C(C[C](F)F)=C(F)F(7490)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {9,S} {10,D}
8  C u1 p0 c0 {1,S} {2,S} {6,S}
9  C u1 p0 c0 {3,S} {7,S} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {7,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-687.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506404,0.0821003,-9.63391e-05,5.93859e-08,-1.49101e-11,-82529.6,31.714], Tmin=(100,'K'), Tmax=(958.197,'K')), NASAPolynomial(coeffs=[12.8859,0.0304206,-1.54351e-05,3.09515e-09,-2.23027e-13,-84901.9,-27.4807], Tmin=(958.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cs-CdHH)(F1s)(F1s)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C]F(744)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = '[CH2]C(F)(F)[C]=CF(1994)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u1 p0 c0 {4,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {7,D} {10,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-172.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,3000,3100,440,815,1455,1000,615,860,1140,1343,3152,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.277176,'amu*angstrom^2'), symmetry=1, barrier=(6.37281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278548,'amu*angstrom^2'), symmetry=1, barrier=(6.40437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19016,0.0690239,-0.000114513,1.04213e-07,-3.7441e-11,-20606.5,23.3252], Tmin=(100,'K'), Tmax=(797.093,'K')), NASAPolynomial(coeffs=[7.46278,0.0265528,-1.39007e-05,2.76061e-09,-1.94585e-13,-21257.2,-3.32326], Tmin=(797.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCd)(H)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = '[C]=CF(1054)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u2 p0 c0 {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (404.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([194,682,905,1196,1383,3221],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62308,0.0070472,-1.17409e-06,-1.98501e-09,8.12281e-13,48709.9,8.54797], Tmin=(100,'K'), Tmax=(1284.59,'K')), NASAPolynomial(coeffs=[5.40185,0.00468,-2.11337e-06,4.24439e-10,-3.06832e-14,47991.3,-1.49755], Tmin=(1284.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'F[C](F)C[C](F)F(795)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 C u1 p0 c0 {3,S} {4,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-548.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,146,234,414,562,504,606,1176,1296,1354,1460,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.12056,'amu*angstrom^2'), symmetry=1, barrier=(2.77191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119179,'amu*angstrom^2'), symmetry=1, barrier=(2.74016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32793,0.0709008,-0.000138892,1.37932e-07,-5.10408e-11,-65872.1,22.2883], Tmin=(100,'K'), Tmax=(850.778,'K')), NASAPolynomial(coeffs=[4.5236,0.0283284,-1.52638e-05,2.99877e-09,-2.07429e-13,-65418.9,13.2467], Tmin=(850.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FC=C1C(F)(F)CC1(F)F(7696)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-930.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896064,0.0764497,-8.71687e-05,5.58957e-08,-1.53246e-11,-111861,21.8337], Tmin=(100,'K'), Tmax=(860.814,'K')), NASAPolynomial(coeffs=[9.14322,0.0381269,-2.03895e-05,4.17754e-09,-3.04339e-13,-113281,-16.7179], Tmin=(860.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-930.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FC=CC(F)(F)C=C(F)F(7697)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,D} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {7,D}
10 C u0 p0 c0 {5,S} {8,D} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-925.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280392,0.0892259,-0.000128339,1.02278e-07,-3.3306e-11,-111227,27.9784], Tmin=(100,'K'), Tmax=(746.795,'K')), NASAPolynomial(coeffs=[10.5625,0.0341543,-1.77263e-05,3.53669e-09,-2.52012e-13,-112763,-18.625], Tmin=(746.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-925.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH)"""),
)

species(
    label = 'F[CH][C]=C(F)F(7343)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-200.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.805],'cm^-1')),
        HinderedRotor(inertia=(0.355657,'amu*angstrom^2'), symmetry=1, barrier=(8.17725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377451,-4.40201e-05,2.68133e-08,-6.58282e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.946,'K')), NASAPolynomial(coeffs=[8.46194,0.0130101,-6.31221e-06,1.26455e-09,-9.14158e-14,-25275.2,-12.1013], Tmin=(983.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'FC=[C]C(F)(F)C=C(F)F(7698)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-688.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2329,'amu*angstrom^2'), symmetry=1, barrier=(28.3469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23147,'amu*angstrom^2'), symmetry=1, barrier=(28.314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0644733,0.0962131,-0.000162631,1.45663e-07,-5.1156e-11,-82615.7,29.2123], Tmin=(100,'K'), Tmax=(809.074,'K')), NASAPolynomial(coeffs=[9.93846,0.0326351,-1.73913e-05,3.43652e-09,-2.40684e-13,-83730.3,-13.3462], Tmin=(809.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-688.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Cds_S)"""),
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
    label = 'FC=C=C(F)C[C](F)F(7699)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-449.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,145,326,398,834,1303,190,488,555,1236,1407,113,247,382,1207,3490,540,610,2055,300.02,785.16,1359.1],'cm^-1')),
        HinderedRotor(inertia=(0.00187248,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0018727,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862589,0.0746733,-9.46363e-05,6.50145e-08,-1.82775e-11,-53981.8,27.7401], Tmin=(100,'K'), Tmax=(859.268,'K')), NASAPolynomial(coeffs=[10.7131,0.0288183,-1.4589e-05,2.90978e-09,-2.08503e-13,-55674.6,-18.2888], Tmin=(859.268,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-449.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CdCddCF) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2][C](F)F(507)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25818,0.018626,-1.56415e-05,8.22322e-09,-2.04586e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.641,'K')), NASAPolynomial(coeffs=[4.47008,0.0131648,-6.41281e-06,1.29207e-09,-9.37479e-14,-13305.9,7.85294], Tmin=(887.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'C#CC(F)(F)C[C](F)F(7700)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u0 p0 c0 {6,S} {9,T}
9  C u0 p0 c0 {8,T} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-511.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,154,355,414,641,686,1150,1196,190,488,555,1236,1407,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.291055,'amu*angstrom^2'), symmetry=1, barrier=(6.69192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.40271,'amu*angstrom^2'), symmetry=1, barrier=(55.243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5758,'amu*angstrom^2'), symmetry=1, barrier=(13.2388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.108614,0.0941582,-0.000156325,1.35689e-07,-4.58025e-11,-61415.3,26.7638], Tmin=(100,'K'), Tmax=(840.54,'K')), NASAPolynomial(coeffs=[10.7362,0.0291401,-1.45215e-05,2.77559e-09,-1.90164e-13,-62691.7,-19.6271], Tmin=(840.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-511.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FC#CC(F)(F)C[C](F)F(7701)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  C u0 p0 c0 {7,S} {10,T}
10 C u0 p0 c0 {5,S} {9,T}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-611.498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,154,355,414,641,686,1150,1196,190,488,555,1236,1407,2175,525,239,401,1367,180,180,2510.58],'cm^-1')),
        HinderedRotor(inertia=(0.274065,'amu*angstrom^2'), symmetry=1, barrier=(6.30129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27198,'amu*angstrom^2'), symmetry=1, barrier=(6.25336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92529,'amu*angstrom^2'), symmetry=1, barrier=(44.2661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.117605,0.10324,-0.000186147,1.71892e-07,-6.0645e-11,-73410.2,29.8958], Tmin=(100,'K'), Tmax=(844.877,'K')), NASAPolynomial(coeffs=[9.0924,0.0341095,-1.80919e-05,3.52343e-09,-2.42974e-13,-74055.4,-7.59292], Tmin=(844.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Ct-CtCs) + group(CtCF) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'FC=[C]C(F)=C[C](F)F(7702)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,D} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {5,D} {9,S}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-336.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,86,203,488,582,605,741,161,297,490,584,780,1358,615,860,1140,1343,3152,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.00425008,'amu*angstrom^2'), symmetry=1, barrier=(4.40065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30238,'amu*angstrom^2'), symmetry=1, barrier=(52.9363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.872645,0.0751001,-0.000103232,7.74178e-08,-2.37407e-11,-40350.6,27.47], Tmin=(100,'K'), Tmax=(790.772,'K')), NASAPolynomial(coeffs=[10.0705,0.0285728,-1.49725e-05,3.0081e-09,-2.1569e-13,-41805.3,-14.745], Tmin=(790.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cd-F1sCd)(Cd-F1sH))"""),
)

species(
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'FC=[C]C(F)(F)[CH]C(F)F(7703)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-616.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0502715,0.0992243,-0.000165005,1.49654e-07,-5.38876e-11,-73980.5,33.3081], Tmin=(100,'K'), Tmax=(785.383,'K')), NASAPolynomial(coeffs=[9.29456,0.0375632,-2.03717e-05,4.0772e-09,-2.88519e-13,-75014.5,-6.75564], Tmin=(785.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-616.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cs_S) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'F[C]=CC(F)(F)C[C](F)F(7704)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {13,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-619.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.268406,0.105144,-0.000180904,1.64097e-07,-5.79979e-11,-74348.5,31.8428], Tmin=(100,'K'), Tmax=(817.53,'K')), NASAPolynomial(coeffs=[9.94842,0.0361075,-1.92879e-05,3.80621e-09,-2.65991e-13,-75382.5,-11.4961], Tmin=(817.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC=CC(F)(F)[CH][C](F)F(7705)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u1 p0 c0 {6,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u0 p0 c0 {5,S} {8,D} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-674.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0492056,0.100984,-0.000175167,1.63293e-07,-5.93353e-11,-81041.3,33.1198], Tmin=(100,'K'), Tmax=(811.621,'K')), NASAPolynomial(coeffs=[8.14878,0.0391171,-2.11575e-05,4.20483e-09,-2.95327e-13,-81665.1,-0.364424], Tmin=(811.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-674.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]=[C]C(F)(F)CC(F)F(7706)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {13,S}
9  C u1 p0 c0 {7,S} {10,D}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-560.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,136,307,446,511,682,757,1180,1185,235,523,627,1123,1142,1372,1406,3097,1685,370,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.265836,'amu*angstrom^2'), symmetry=1, barrier=(6.11209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265844,'amu*angstrom^2'), symmetry=1, barrier=(6.11227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.620165,'amu*angstrom^2'), symmetry=1, barrier=(14.2588,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.266652,0.103344,-0.000170564,1.50152e-07,-5.23726e-11,-67287.9,32.0214], Tmin=(100,'K'), Tmax=(788.782,'K')), NASAPolynomial(coeffs=[11.097,0.0345488,-1.84995e-05,3.67793e-09,-2.5913e-13,-68733.1,-17.9028], Tmin=(788.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[CH][C]=C(F)CC(F)(F)F(7707)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u1 p0 c0 {5,S} {10,S} {13,S}
10 C u1 p0 c0 {8,D} {9,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-706.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442383,0.083353,-9.9521e-05,6.19711e-08,-1.56359e-11,-84879.9,29.8129], Tmin=(100,'K'), Tmax=(956.158,'K')), NASAPolynomial(coeffs=[13.3307,0.0294355,-1.49361e-05,2.99524e-09,-2.15831e-13,-87344.6,-31.788], Tmin=(956.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-706.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-CsF1s))"""),
)

species(
    label = 'FC=C(F)[C](F)C[C](F)F(7708)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
7  C u1 p0 c0 {1,S} {6,S} {8,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 C u0 p0 c0 {5,S} {8,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-677.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.215398,0.101822,-0.000165454,1.44417e-07,-4.99932e-11,-81321.1,31.6139], Tmin=(100,'K'), Tmax=(794.784,'K')), NASAPolynomial(coeffs=[10.9459,0.0347192,-1.81838e-05,3.58312e-09,-2.51237e-13,-82750,-17.4974], Tmin=(794.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-677.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCdCsF1s) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH]=[C]C(F)(F)CC(F)(F)F(7709)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {10,D}
10 C u1 p0 c0 {9,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-635.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0411751,0.0965944,-0.000141429,1.04559e-07,-2.82847e-11,-76311.9,30.3365], Tmin=(100,'K'), Tmax=(628.419,'K')), NASAPolynomial(coeffs=[12.0138,0.0324799,-1.72581e-05,3.45451e-09,-2.45817e-13,-78055.5,-23.7629], Tmin=(628.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-635.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFF) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sCs)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)C[C](F)F(7710)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u1 p0 c0 {8,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-611.285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,2750,2850,1437.5,1250,1305,750,350,246,474,533,1155,190,488,555,1236,1407,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.149501,'amu*angstrom^2'), symmetry=1, barrier=(3.43733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152642,'amu*angstrom^2'), symmetry=1, barrier=(3.50954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150911,'amu*angstrom^2'), symmetry=1, barrier=(3.46975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.387909,0.107301,-0.00018326,1.63581e-07,-5.68892e-11,-73372.9,31.681], Tmin=(100,'K'), Tmax=(820.993,'K')), NASAPolynomial(coeffs=[10.9543,0.0345126,-1.82476e-05,3.58195e-09,-2.4945e-13,-74644.6,-17.2038], Tmin=(820.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-CdsHH) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    E0 = (-247.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-81.1574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-122.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (223.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (218.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-239.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-158.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-165.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-82.1002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (35.4021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-99.4945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-9.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-37.4204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (65.9504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-27.0105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-11.5575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-67.4686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-43.8304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-105.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-165.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-66.4306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-97.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-62.2837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-66.5501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['CH2CF2(57)', 'FC=C=C(F)F(5948)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C(F)(F)[C]=CF(7489)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['F[CH]C(C[C](F)F)=C(F)F(7490)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]F(744)', '[CH2]C(F)(F)[C]=CF(1994)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=CF(1054)', 'F[C](F)C[C](F)F(795)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['FC=C1C(F)(F)CC1(F)F(7696)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['FC=CC(F)(F)C=C(F)F(7697)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CF2(57)', 'F[CH][C]=C(F)F(7343)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.000388714,'m^3/(mol*s)'), n=2.45366, Ea=(34.6812,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.450791267159554, var=0.5685704571032071, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Sp-5R!H=4R!H_Ext-2R!H-R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Sp-5R!H=4R!H_Ext-2R!H-R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'FC=[C]C(F)(F)C=C(F)F(7698)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(306600,'m^3/(mol*s)'), n=0.481, Ea=(31.8287,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'FC=C=C(F)C[C](F)F(7699)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(49.9654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C](F)F(507)', 'FC=C=C(F)F(5948)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(11.9821,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'C#CC(F)(F)C[C](F)F(7700)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(67.356,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'FC#CC(F)(F)C[C](F)F(7701)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](F)F(507)', 'F[CH][C]=C(F)F(7343)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(13.391,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'FC=[C]C(F)=C[C](F)F(7702)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(228.213,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2(43)', '[CH2]C(F)(F)[C]=CF(1994)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(1.99695,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['FC=[C]C(F)(F)[CH]C(F)F(7703)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.76836e+09,'s^-1'), n=1.1815, Ea=(180.499,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['F[C]=CC(F)(F)C[C](F)F(7704)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['FC=CC(F)(F)[CH][C](F)F(7705)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C]=[C]C(F)(F)CC(F)F(7706)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['F[CH][C]=C(F)CC(F)(F)F(7707)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(181.537,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['FC=C(F)[C](F)C[C](F)F(7708)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(150.509,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    products = ['[CH]=[C]C(F)(F)CC(F)(F)F(7709)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(185.684,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(F)C(F)(F)C[C](F)F(7710)'],
    products = ['FC=[C]C(F)(F)C[C](F)F(7492)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(182.461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #2749',
    isomers = [
        'FC=[C]C(F)(F)C[C](F)F(7492)',
    ],
    reactants = [
        ('CH2CF2(57)', 'FC=C=C(F)F(5948)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2749',
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

