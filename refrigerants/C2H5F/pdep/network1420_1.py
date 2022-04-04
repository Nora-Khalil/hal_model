species(
    label = 'CCOOC=C(F)O[O](3746)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {3,S} {9,D} {15,S}
9  C u0 p0 c0 {1,S} {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-170.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.285407,'amu*angstrom^2'), symmetry=1, barrier=(6.56208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285661,'amu*angstrom^2'), symmetry=1, barrier=(6.5679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286137,'amu*angstrom^2'), symmetry=1, barrier=(6.57884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285577,'amu*angstrom^2'), symmetry=1, barrier=(6.56597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125387,0.0961713,-0.000153365,1.42291e-07,-5.24254e-11,-20374.8,33.9277], Tmin=(100,'K'), Tmax=(806.533,'K')), NASAPolynomial(coeffs=[6.2125,0.0461944,-2.36156e-05,4.62313e-09,-3.23324e-13,-20713.1,9.85961], Tmin=(806.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)C=O(1780)',
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
        HinderedRotor(inertia=(1.18318,'amu*angstrom^2'), symmetry=1, barrier=(27.2035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603053,'amu*angstrom^2'), symmetry=1, barrier=(13.8654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3728.91,'J/mol'), sigma=(5.84942,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.45 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15437,0.0440803,-5.80203e-05,4.31468e-08,-1.32604e-11,-39422,18.3048], Tmin=(100,'K'), Tmax=(786.661,'K')), NASAPolynomial(coeffs=[7.18331,0.0185099,-9.26421e-06,1.8289e-09,-1.29954e-13,-40213.3,-4.75028], Tmin=(786.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.65,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'COOC=C(F)O[O](3464)',
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
        HinderedRotor(inertia=(0.313215,'amu*angstrom^2'), symmetry=1, barrier=(7.20143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.313323,'amu*angstrom^2'), symmetry=1, barrier=(7.20392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.313545,'amu*angstrom^2'), symmetry=1, barrier=(7.20901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.313488,'amu*angstrom^2'), symmetry=1, barrier=(7.20772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.06,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3878.7,'J/mol'), sigma=(6.18343,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.84 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916041,0.0794773,-0.000138685,1.35829e-07,-5.12082e-11,-15888,29.2588], Tmin=(100,'K'), Tmax=(826.685,'K')), NASAPolynomial(coeffs=[4.0038,0.0395014,-2.07236e-05,4.06761e-09,-2.83569e-13,-15543,20.124], Tmin=(826.685,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = 'CCOC(F)(C=O)O[O](4301)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {9,D}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  C u0 p0 c0 {4,D} {7,S} {15,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-553.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.219689,0.101572,-0.000156569,1.35062e-07,-4.66455e-11,-66433.8,30.7287], Tmin=(100,'K'), Tmax=(800.159,'K')), NASAPolynomial(coeffs=[9.98214,0.0395616,-1.9682e-05,3.81396e-09,-2.65489e-13,-67713.9,-14.0122], Tmin=(800.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-553.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CCOO[CH]C(=O)F(4302)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  C u1 p0 c0 {3,S} {8,S} {14,S}
8  C u0 p0 c0 {1,S} {4,D} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-404.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,180,3426.74],'cm^-1')),
        HinderedRotor(inertia=(1.80493,'amu*angstrom^2'), symmetry=1, barrier=(41.4988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.612525,'amu*angstrom^2'), symmetry=1, barrier=(14.0832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80443,'amu*angstrom^2'), symmetry=1, barrier=(41.4873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80493,'amu*angstrom^2'), symmetry=1, barrier=(41.4989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80441,'amu*angstrom^2'), symmetry=1, barrier=(41.4869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882473,0.0750053,-8.1009e-05,5.0264e-08,-1.33484e-11,-48577.8,26.7704], Tmin=(100,'K'), Tmax=(888.498,'K')), NASAPolynomial(coeffs=[9.03178,0.0383171,-1.90701e-05,3.7892e-09,-2.71525e-13,-50025.9,-11.5818], Tmin=(888.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(OCJC=O)"""),
)

species(
    label = 'FC1=COOO1(3487)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45638,-0.00295965,7.42287e-05,-9.36153e-08,3.49541e-11,-7551.03,17.4586], Tmin=(100,'K'), Tmax=(965.702,'K')), NASAPolynomial(coeffs=[9.89832,0.00855837,-2.9985e-06,6.60683e-10,-5.63558e-14,-10576.5,-22.6178], Tmin=(965.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.0576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclopentane)"""),
)

species(
    label = 'CC[O](149)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-25.0859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,312.548],'cm^-1')),
        HinderedRotor(inertia=(0.118919,'amu*angstrom^2'), symmetry=1, barrier=(2.73417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0605,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26906,0.00933563,2.96317e-05,-4.53411e-08,1.88796e-11,-2950.23,10.4201], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.55054,0.0132526,-4.74726e-06,7.64699e-10,-4.57008e-14,-4471.92,-9.61231], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-25.0859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH3CH2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CCOOC1OO[C]1F(4303)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
7  C u0 p0 c0 {2,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
9  C u1 p0 c0 {1,S} {5,S} {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (-210.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.804032,0.0778268,-8.16996e-05,4.99056e-08,-1.32613e-11,-25184.7,30.3975], Tmin=(100,'K'), Tmax=(880.164,'K')), NASAPolynomial(coeffs=[8.57652,0.042504,-2.15019e-05,4.31002e-09,-3.10471e-13,-26552.9,-6.1082], Tmin=(880.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'CCOO[CH]C1(F)OO1(4304)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
7  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
9  C u1 p0 c0 {5,S} {6,S} {15,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-251.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462622,0.0849591,-9.63766e-05,6.2515e-08,-1.71881e-11,-30130.8,29.1904], Tmin=(100,'K'), Tmax=(864.404,'K')), NASAPolynomial(coeffs=[9.89387,0.0413154,-2.06398e-05,4.10211e-09,-2.93789e-13,-31761.2,-14.9354], Tmin=(864.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(O2s-O2s-Cs(F)(C)) + radical(CCsJOOC)"""),
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
    label = 'CCOOC=[C]O[O](4305)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (227.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.366537,'amu*angstrom^2'), symmetry=1, barrier=(8.42741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366477,'amu*angstrom^2'), symmetry=1, barrier=(8.42603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366575,'amu*angstrom^2'), symmetry=1, barrier=(8.42829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36649,'amu*angstrom^2'), symmetry=1, barrier=(8.42634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11683,'amu*angstrom^2'), symmetry=1, barrier=(25.678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.529251,0.0817383,-0.000110059,8.44135e-08,-2.65024e-11,27503.4,33.6029], Tmin=(100,'K'), Tmax=(774.46,'K')), NASAPolynomial(coeffs=[9.99043,0.0328735,-1.54182e-05,2.9469e-09,-2.05054e-13,26037.9,-9.62388], Tmin=(774.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]O[C](F)C=O(1865)',
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
        HinderedRotor(inertia=(4.19309,'amu*angstrom^2'), symmetry=1, barrier=(96.4073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.19306,'amu*angstrom^2'), symmetry=1, barrier=(96.4068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35164,0.0412356,-6.68966e-05,6.24094e-08,-2.28577e-11,-22685.8,17.8038], Tmin=(100,'K'), Tmax=(824.81,'K')), NASAPolynomial(coeffs=[4.7315,0.019829,-1.00256e-05,1.94133e-09,-1.34579e-13,-22742.8,8.81502], Tmin=(824.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]OC=C(F)O[O](3491)',
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
        HinderedRotor(inertia=(0.234516,'amu*angstrom^2'), symmetry=1, barrier=(5.39198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68878,0.0630762,-0.000129041,1.29579e-07,-4.76621e-11,2374.86,24.9913], Tmin=(100,'K'), Tmax=(871.39,'K')), NASAPolynomial(coeffs=[3.58343,0.0251327,-1.33818e-05,2.57537e-09,-1.74839e-13,3155.04,22.4828], Tmin=(871.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.1511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = 'C2H5(32)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.89,1642.9,3622.3,3622.32],'cm^-1')),
        HinderedRotor(inertia=(0.866818,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=C(F)O[O](1331)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00596,0.051545,-0.000104345,9.71569e-08,-3.31602e-11,19053.1,15.8311], Tmin=(100,'K'), Tmax=(897.32,'K')), NASAPolynomial(coeffs=[7.11387,0.0105561,-5.37029e-06,9.95723e-10,-6.48099e-14,18869.9,-4.17112], Tmin=(897.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'CCO[O](193)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-35.8481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.196253,'amu*angstrom^2'), symmetry=1, barrier=(4.51224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197956,'amu*angstrom^2'), symmetry=1, barrier=(4.5514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88894,0.0185233,1.04542e-05,-2.35199e-08,9.42774e-12,-4266.36,15.0609], Tmin=(100,'K'), Tmax=(1009.3,'K')), NASAPolynomial(coeffs=[7.15149,0.0156281,-6.0467e-06,1.12072e-09,-7.93638e-14,-5839.77,-9.07494], Tmin=(1009.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.8481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CCOOC=[C]F(4306)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {7,D} {13,S}
7  C u1 p0 c0 {1,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (14.1933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.296857,'amu*angstrom^2'), symmetry=1, barrier=(6.82532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7476,'amu*angstrom^2'), symmetry=1, barrier=(17.1888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297211,'amu*angstrom^2'), symmetry=1, barrier=(6.83345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.307,'amu*angstrom^2'), symmetry=1, barrier=(53.0424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3922,0.0658144,-6.56013e-05,1.85383e-08,1.32554e-11,1792.73,24.6147], Tmin=(100,'K'), Tmax=(539.285,'K')), NASAPolynomial(coeffs=[7.18144,0.0357268,-1.76628e-05,3.46928e-09,-2.45717e-13,981.43,-1.47288], Tmin=(539.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.1933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    label = '[CH2]OOC=C(F)O[O](3493)',
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
        HinderedRotor(inertia=(0.228242,'amu*angstrom^2'), symmetry=1, barrier=(5.24774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227925,'amu*angstrom^2'), symmetry=1, barrier=(5.24044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228284,'amu*angstrom^2'), symmetry=1, barrier=(5.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.39235,'amu*angstrom^2'), symmetry=1, barrier=(55.0049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864033,0.0797895,-0.000140546,1.34655e-07,-4.97472e-11,7399.04,29.1424], Tmin=(100,'K'), Tmax=(825.995,'K')), NASAPolynomial(coeffs=[5.57945,0.0347714,-1.85092e-05,3.64432e-09,-2.54233e-13,7376.79,11.8755], Tmin=(825.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.6674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(CsJOOC)"""),
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
    label = 'C[CH]OOC=C(F)O[O](4307)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {2,S} {6,S} {13,S}
8  C u0 p0 c0 {3,S} {9,D} {14,S}
9  C u0 p0 c0 {1,S} {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (15.9635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.165479,'amu*angstrom^2'), symmetry=1, barrier=(3.80468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165638,'amu*angstrom^2'), symmetry=1, barrier=(3.80834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164746,'amu*angstrom^2'), symmetry=1, barrier=(3.78782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126767,0.0980891,-0.000167712,1.59638e-07,-5.89161e-11,2046.93,34.2176], Tmin=(100,'K'), Tmax=(826.137,'K')), NASAPolynomial(coeffs=[5.67634,0.0448964,-2.33378e-05,4.56402e-09,-3.17609e-13,2028.26,13.9406], Tmin=(826.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.9635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]COOC=C(F)O[O](4308)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {9,D} {12,S}
8  C u1 p0 c0 {6,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {4,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (43.4853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.223702,'amu*angstrom^2'), symmetry=1, barrier=(5.14335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223924,'amu*angstrom^2'), symmetry=1, barrier=(5.14846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223696,'amu*angstrom^2'), symmetry=1, barrier=(5.14322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223757,'amu*angstrom^2'), symmetry=1, barrier=(5.1446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.1838,'amu*angstrom^2'), symmetry=1, barrier=(73.2018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0285459,0.0992517,-0.000167692,1.5613e-07,-5.66649e-11,5361.57,36.2708], Tmin=(100,'K'), Tmax=(822.88,'K')), NASAPolynomial(coeffs=[7.11244,0.0420175,-2.18008e-05,4.26366e-09,-2.96846e-13,4967.64,8.16635], Tmin=(822.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.4853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'CCOO[C]=C(F)O[O](4309)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (69.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,293,496,537,1218,1685,370,180,1466.2],'cm^-1')),
        HinderedRotor(inertia=(0.242877,'amu*angstrom^2'), symmetry=1, barrier=(5.58423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242889,'amu*angstrom^2'), symmetry=1, barrier=(5.5845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242897,'amu*angstrom^2'), symmetry=1, barrier=(5.58467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243012,'amu*angstrom^2'), symmetry=1, barrier=(5.58732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242877,'amu*angstrom^2'), symmetry=1, barrier=(5.58423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.297133,0.0983936,-0.000181907,1.83183e-07,-6.94654e-11,8447.83,36.4385], Tmin=(100,'K'), Tmax=(842.817,'K')), NASAPolynomial(coeffs=[2.28329,0.0499481,-2.62422e-05,5.1218e-09,-3.54585e-13,9498.88,35.4176], Tmin=(842.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
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
    label = 'CCOOC#CO[O](4310)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,T}
8  C u0 p0 c0 {3,S} {7,T}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (278.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,180],'cm^-1')),
        HinderedRotor(inertia=(0.71065,'amu*angstrom^2'), symmetry=1, barrier=(16.3392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710653,'amu*angstrom^2'), symmetry=1, barrier=(16.3393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04891,'amu*angstrom^2'), symmetry=1, barrier=(47.1085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08208,'amu*angstrom^2'), symmetry=1, barrier=(24.8793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04649,'amu*angstrom^2'), symmetry=1, barrier=(47.0528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544923,0.0824666,-0.000114937,9.12692e-08,-2.98497e-11,33657.8,29.6341], Tmin=(100,'K'), Tmax=(742.169,'K')), NASAPolynomial(coeffs=[9.53221,0.0340325,-1.70545e-05,3.35125e-09,-2.36808e-13,32323.7,-11.0452], Tmin=(742.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'CCOO[C]=C(F)OO(4311)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {8,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u0 p2 c0 {3,S} {15,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {1,S} {3,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-82.7377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,293,496,537,1218,1685,370,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0462078,0.101853,-0.000176903,1.73891e-07,-6.59898e-11,-9823.15,36.5794], Tmin=(100,'K'), Tmax=(821.581,'K')), NASAPolynomial(coeffs=[3.81336,0.0517103,-2.7293e-05,5.37682e-09,-3.75911e-13,-9368.85,25.6773], Tmin=(821.581,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.7377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]OOC=C(F)OO(4312)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {15,S}
6  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {2,S} {6,S} {13,S}
8  C u0 p0 c0 {3,S} {9,D} {14,S}
9  C u0 p0 c0 {1,S} {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-136.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.104681,0.101309,-0.000161845,1.49244e-07,-5.50138e-11,-16224.9,34.2893], Tmin=(100,'K'), Tmax=(791.685,'K')), NASAPolynomial(coeffs=[7.11534,0.0468211,-2.44856e-05,4.84256e-09,-3.40922e-13,-16803.7,4.70778], Tmin=(791.685,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]COOC=C(F)OO(4313)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {15,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {8,D} {14,S}
8  C u0 p0 c0 {1,S} {4,S} {7,D}
9  C u1 p0 c0 {6,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
"""),
    E0 = (-108.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,326,540,652,719,1357,3000,3100,440,815,1455,1000,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979728,0.0826291,-5.01454e-05,-1.083e-07,1.47825e-10,-12960.2,32.4048], Tmin=(100,'K'), Tmax=(449.6,'K')), NASAPolynomial(coeffs=[8.89245,0.0432936,-2.2544e-05,4.44106e-09,-3.11412e-13,-13985.6,-2.93556], Tmin=(449.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CJCOOH)"""),
)

species(
    label = 'CCOOC=[C]OOF(4314)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {1,S} {4,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {3,S} {9,D} {15,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
"""),
    E0 = (170.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.886871,'amu*angstrom^2'), symmetry=1, barrier=(20.3909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.884285,'amu*angstrom^2'), symmetry=1, barrier=(20.3315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890333,'amu*angstrom^2'), symmetry=1, barrier=(20.4705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.883239,'amu*angstrom^2'), symmetry=1, barrier=(20.3074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0937593,0.0925905,-0.000118765,8.38036e-08,-2.42129e-11,20698.3,35.8008], Tmin=(100,'K'), Tmax=(837.512,'K')), NASAPolynomial(coeffs=[11.8851,0.0362765,-1.791e-05,3.52494e-09,-2.50434e-13,18723.1,-18.9954], Tmin=(837.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
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
    E0 = (-141.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (347.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-83.5395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-159.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-86.2278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-42.3259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-130.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (302.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-160.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (128.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (123.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (7.48215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (198.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (229.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (257.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (282.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (285.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-36.5134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-92.5382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (34.3972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (265.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCOOC=C(F)O[O](3746)'],
    products = ['[O]OC(F)C=O(1780)', 'CH3CHO(36)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(702.966,'s^-1'), n=2.72887, Ea=(27.2136,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', 'COOC=C(F)O[O](3464)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.26474e-05,'m^3/(mol*s)'), n=2.95311, Ea=(59.8418,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CCOOC=C(F)O[O](3746)'],
    products = ['CCOC(F)(C=O)O[O](4301)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(85.022,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'CCOO[CH]C(=O)F(4302)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CCOOC=C(F)O[O](3746)'],
    products = ['FC1=COOO1(3487)', 'CC[O](149)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(82.3337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 80.0 to 82.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CCOOC=C(F)O[O](3746)'],
    products = ['CCOOC1OO[C]1F(4303)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CCOOC=C(F)O[O](3746)'],
    products = ['CCOO[CH]C1(F)OO1(4304)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(38.4928,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'CCOOC=[C]O[O](4305)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]O[C](F)C=O(1865)', 'CC[O](149)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(51.7497,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC=C(F)O[O](3491)', 'C2H5(32)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.01905e+07,'m^3/(mol*s)'), n=-0.198371, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.054859309462099146, var=1.178102591472186, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(F)O[O](1331)', 'CCO[O](193)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O2(2)', 'CCOOC=[C]F(4306)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(19)', '[CH2]OOC=C(F)O[O](3493)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'C[CH]OOC=C(F)O[O](4307)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[CH2]COOC=C(F)O[O](4308)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'CCOO[C]=C(F)O[O](4309)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'CCOOC#CO[O](4310)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(286.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCOO[C]=C(F)OO(4311)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CCOOC=C(F)O[O](3746)'],
    products = ['C[CH]OOC=C(F)OO(4312)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.74e+06,'s^-1'), n=0.99, Ea=(76.0233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 267 used for R7H_OOCCCC(Cs/Cs);O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R7H_OOCCCC(Cs/Cs);O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]COOC=C(F)OO(4313)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.27e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H;C_rad_out_2H;XH_out] for rate rule [R8H;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCOOC=[C]OOF(4314)'],
    products = ['CCOOC=C(F)O[O](3746)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.8166,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1420',
    isomers = [
        'CCOOC=C(F)O[O](3746)',
    ],
    reactants = [
        ('[O]OC(F)C=O(1780)', 'CH3CHO(36)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1420',
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

