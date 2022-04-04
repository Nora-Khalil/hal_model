species(
    label = '[CH2]OO[CH]C(=O)F(6174)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u1 p0 c0 {2,S} {6,S} {8,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u1 p0 c0 {3,S} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-173.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,611,648,830,1210,1753,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.597675,'amu*angstrom^2'), symmetry=1, barrier=(45.6269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98455,'amu*angstrom^2'), symmetry=1, barrier=(45.6286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98448,'amu*angstrom^2'), symmetry=1, barrier=(45.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.01759,'amu*angstrom^2'), symmetry=1, barrier=(45.6276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51557,0.0599487,-7.32694e-05,5.00336e-08,-1.42964e-11,-20799.5,22.3588], Tmin=(100,'K'), Tmax=(837.815,'K')), NASAPolynomial(coeffs=[8.46047,0.0267917,-1.39061e-05,2.7971e-09,-2.01353e-13,-21963.2,-9.91742], Tmin=(837.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsHHH) + group(COCsFO) + radical(OCJC=O) + radical(CsJOOC)"""),
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
    label = 'O=CC(=O)F(335)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]O[O](231)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (206.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,291.401,3698.4],'cm^-1')),
        HinderedRotor(inertia=(0.37416,'amu*angstrom^2'), symmetry=1, barrier=(22.5475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92061,0.0189601,-2.08419e-05,1.10972e-08,-2.15337e-12,24926,11.4734], Tmin=(100,'K'), Tmax=(1506.17,'K')), NASAPolynomial(coeffs=[8.0388,0.00135553,6.86083e-07,-2.00091e-10,1.53442e-14,23839.3,-13.8045], Tmin=(1506.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = '[CH]C(=O)F(323)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-5.0725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([486,617,768,1157,1926,180,1655.05,1655.52],'cm^-1')),
        HinderedRotor(inertia=(0.0191536,'amu*angstrom^2'), symmetry=1, barrier=(5.31745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32766,0.0152939,-1.10758e-05,3.91582e-09,-5.6142e-13,-586.537,12.101], Tmin=(100,'K'), Tmax=(1580.41,'K')), NASAPolynomial(coeffs=[6.53552,0.00717479,-3.36979e-06,6.65147e-10,-4.72043e-14,-1600.48,-4.84321], Tmin=(1580.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]O[CH]C(=O)F(1948)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u1 p0 c0 {2,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-215.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(0.2824,'amu*angstrom^2'), symmetry=1, barrier=(6.49293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4819,'amu*angstrom^2'), symmetry=1, barrier=(57.0637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06408,0.0468221,-7.59984e-05,6.61717e-08,-2.2668e-11,-25811.9,19.1781], Tmin=(100,'K'), Tmax=(821.51,'K')), NASAPolynomial(coeffs=[7.05431,0.0160915,-8.14163e-06,1.57328e-09,-1.08835e-13,-26414.8,-2.59461], Tmin=(821.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(ROOJ) + radical(OCJC=O)"""),
)

species(
    label = 'O=C(F)C1COO1(6176)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-386.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9661,0.0380115,-1.91932e-05,-6.93064e-11,1.86995e-12,-46386.9,19.7015], Tmin=(100,'K'), Tmax=(1191.33,'K')), NASAPolynomial(coeffs=[11.1326,0.0182815,-8.26102e-06,1.59701e-09,-1.13195e-13,-49355,-29.4167], Tmin=(1191.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-386.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(COCsFO) + ring(12dioxetane)"""),
)

species(
    label = '[CH2]OOC1O[C]1F(6639)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
6  C u1 p0 c0 {1,S} {2,S} {5,S}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (18.4157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6608,0.0500633,-4.17876e-05,1.69359e-08,-2.76967e-12,2299.94,21.0021], Tmin=(100,'K'), Tmax=(1428,'K')), NASAPolynomial(coeffs=[12.6244,0.0193526,-9.52809e-06,1.87522e-09,-1.32962e-13,-831.208,-35.7966], Tmin=(1428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.4157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(Cs-OsHHH) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring) + radical(CsJOOC)"""),
)

species(
    label = 'F[C]1[CH]OOCO1(6640)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
6  C u1 p0 c0 {1,S} {2,S} {7,S}
7  C u1 p0 c0 {4,S} {6,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-126.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39507,0.031068,6.52317e-05,-1.42859e-07,7.00523e-11,-15077.8,18.7613], Tmin=(100,'K'), Tmax=(861.333,'K')), NASAPolynomial(coeffs=[27.6114,-0.0156395,1.58906e-05,-3.52215e-09,2.5259e-13,-22377.6,-119.963], Tmin=(861.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(124trioxane) + radical(CsCsF1sO2s) + radical(CCsJOOC)"""),
)

species(
    label = '[O]C1(F)[CH]OOC1(6641)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {3,S} {5,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-120.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03635,0.0521513,-4.98973e-05,2.48814e-08,-4.67759e-12,-14356.5,21.0868], Tmin=(100,'K'), Tmax=(1536.36,'K')), NASAPolynomial(coeffs=[12.7791,0.0122441,-1.82119e-06,6.53737e-11,4.00974e-15,-16863.1,-37.0223], Tmin=(1536.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-120.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(12dioxolane) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2][O](2599)',
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
    label = '[O][CH]C(=O)F(509)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.101,381.695],'cm^-1')),
        HinderedRotor(inertia=(0.482775,'amu*angstrom^2'), symmetry=1, barrier=(49.7784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80088,0.0290439,-3.02832e-05,1.66615e-08,-3.81631e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.22,'K')), NASAPolynomial(coeffs=[6.95396,0.0128879,-6.71487e-06,1.38086e-09,-1.01081e-13,-26608.7,-6.32755], Tmin=(1028.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
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
    label = '[CH2]OOC=C=O(6642)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u1 p0 c0 {2,S} {8,S} {9,S}
6 C u0 p0 c0 {3,D} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (124.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.72456,'amu*angstrom^2'), symmetry=1, barrier=(16.6591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58423,'amu*angstrom^2'), symmetry=1, barrier=(36.4245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.85273,'amu*angstrom^2'), symmetry=1, barrier=(111.574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82311,0.0482118,-5.20625e-05,2.89278e-08,-6.44161e-12,14998.6,19.6942], Tmin=(100,'K'), Tmax=(1084.98,'K')), NASAPolynomial(coeffs=[10.6123,0.0158089,-7.26529e-06,1.40227e-09,-9.92415e-14,13091.4,-23.4254], Tmin=(1084.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-OsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsJOOC)"""),
)

species(
    label = '[CH2]OOC(F)[C]=O(6630)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u1 p0 c0 {3,S} {9,S} {10,S}
7  C u1 p0 c0 {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-128.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,355,410,600,1181,1341,1420,3056,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.12955,'amu*angstrom^2'), symmetry=1, barrier=(25.9705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40862,'amu*angstrom^2'), symmetry=1, barrier=(9.39499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.39675,'amu*angstrom^2'), symmetry=1, barrier=(55.106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.39649,'amu*angstrom^2'), symmetry=1, barrier=(55.1,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970493,0.0725379,-0.000114897,9.78722e-08,-3.31971e-11,-15320.5,23.4148], Tmin=(100,'K'), Tmax=(791.557,'K')), NASAPolynomial(coeffs=[9.41021,0.024171,-1.24058e-05,2.42555e-09,-1.69479e-13,-16477.5,-14.1976], Tmin=(791.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(CsJOOC) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C][CH]OOCF(6643)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
6  C u1 p0 c0 {3,S} {7,S} {10,S}
7  C u1 p0 c0 {4,D} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-172.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59678,0.0584886,-7.39283e-05,5.42115e-08,-1.68703e-11,-20715.6,23.7284], Tmin=(100,'K'), Tmax=(768.577,'K')), NASAPolynomial(coeffs=[7.40799,0.0282451,-1.49043e-05,3.01464e-09,-2.17453e-13,-21608.9,-2.7778], Tmin=(768.577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFHHO) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]OC([O])(F)C=O(6173)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u0 p0 c0 {4,D} {5,S} {8,S}
7  C u1 p0 c0 {2,S} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-307.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09855,0.0703228,-0.000115614,1.0115e-07,-3.46886e-11,-36935.2,24.216], Tmin=(100,'K'), Tmax=(823.374,'K')), NASAPolynomial(coeffs=[8.66387,0.0234841,-1.19104e-05,2.30639e-09,-1.59664e-13,-37839.1,-8.73585], Tmin=(823.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-307.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsJOCH3)"""),
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
    label = '[CH2]OOC=[C]F(6644)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {6,D} {7,S}
5 C u1 p0 c0 {3,S} {8,S} {9,S}
6 C u1 p0 c0 {1,S} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (245.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(0.177775,'amu*angstrom^2'), symmetry=1, barrier=(4.0874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.09836,'amu*angstrom^2'), symmetry=1, barrier=(94.2294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17994,'amu*angstrom^2'), symmetry=1, barrier=(50.121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67103,0.0563721,-8.6698e-05,7.61217e-08,-2.70258e-11,29586.2,21.3964], Tmin=(100,'K'), Tmax=(778.82,'K')), NASAPolynomial(coeffs=[6.85446,0.0237183,-1.21898e-05,2.39865e-09,-1.68672e-13,28961.8,-1.14055], Tmin=(778.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsJOOC) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'FC1=COOCO1(6645)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-323.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76829,0.0234372,6.91814e-05,-1.29126e-07,5.8062e-11,-38790.4,17.1931], Tmin=(100,'K'), Tmax=(905.75,'K')), NASAPolynomial(coeffs=[24.0774,-0.00739252,8.13362e-06,-1.67913e-09,1.10016e-13,-45608.4,-103.555], Tmin=(905.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsOsHH) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
)

species(
    label = 'FC1=COO1(285)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-52.3112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,326,540,652,719,1357,353.677,1138.56,1138.58,1138.63,1726.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95335,0.00258575,6.07332e-05,-1.08656e-07,5.80728e-11,-6290.59,8.5661], Tmin=(10,'K'), Tmax=(614.296,'K')), NASAPolynomial(coeffs=[2.56312,0.0230016,-1.68658e-05,5.67097e-09,-7.10066e-13,-6334.19,12.8506], Tmin=(614.296,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-52.3112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC1DCOO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O][C](F)C1COO1(6646)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
7  C u1 p0 c0 {1,S} {4,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-16.5431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55547,0.0603498,-9.52817e-05,8.64533e-08,-3.09303e-11,-1908.01,20.3313], Tmin=(100,'K'), Tmax=(828.943,'K')), NASAPolynomial(coeffs=[5.76414,0.027602,-1.35144e-05,2.59046e-09,-1.78695e-13,-2178.39,3.39424], Tmin=(828.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.5431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCFHO) + ring(12dioxetane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = '[CH2]OO[C]=C=O(6635)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u1 p0 c0 {1,S} {7,S} {8,S}
5 C u1 p0 c0 {2,S} {6,D}
6 C u0 p0 c0 {3,D} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (363.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.52922,'amu*angstrom^2'), symmetry=1, barrier=(12.1678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0122591,'amu*angstrom^2'), symmetry=1, barrier=(44.3881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92155,'amu*angstrom^2'), symmetry=1, barrier=(44.1802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90219,0.051592,-8.51021e-05,7.66502e-08,-2.70306e-11,43825.3,22.5343], Tmin=(100,'K'), Tmax=(826.955,'K')), NASAPolynomial(coeffs=[6.40244,0.020057,-1.01849e-05,1.9716e-09,-1.36458e-13,43415,3.69753], Tmin=(826.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-OsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsJOOC) + radical(C=CJO)"""),
)

species(
    label = '[CH2]OO[C]=C(O)F(6647)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u1 p0 c0 {2,S} {8,S} {9,S}
7  C u1 p0 c0 {3,S} {5,D}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (23.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,293,496,537,1218,3000,3100,440,815,1455,1000,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.00447569,'amu*angstrom^2'), symmetry=1, barrier=(5.15839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224243,'amu*angstrom^2'), symmetry=1, barrier=(5.15579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06305,'amu*angstrom^2'), symmetry=1, barrier=(47.4335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06432,'amu*angstrom^2'), symmetry=1, barrier=(47.4628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14522,0.0692441,-0.000111271,9.82185e-08,-3.45723e-11,2972.53,26.4629], Tmin=(100,'K'), Tmax=(789.955,'K')), NASAPolynomial(coeffs=[8.07018,0.0259957,-1.36103e-05,2.68637e-09,-1.88746e-13,2133.78,-3.69715], Tmin=(789.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsJOOC) + radical(C=CJO)"""),
)

species(
    label = 'COO[C]C(=O)F(6648)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,D} {7,S}
7  C u2 p0 c0 {3,S} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-44.7352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926,235.693,235.722],'cm^-1')),
        HinderedRotor(inertia=(0.78606,'amu*angstrom^2'), symmetry=1, barrier=(30.958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785593,'amu*angstrom^2'), symmetry=1, barrier=(30.959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429077,'amu*angstrom^2'), symmetry=1, barrier=(16.8752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17801,'amu*angstrom^2'), symmetry=1, barrier=(7.01768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42816,0.0619848,-7.78159e-05,5.16747e-08,-1.40618e-11,-5292.42,21.9521], Tmin=(100,'K'), Tmax=(884.308,'K')), NASAPolynomial(coeffs=[9.91098,0.0236139,-1.27289e-05,2.60597e-09,-1.89534e-13,-6792.68,-17.9295], Tmin=(884.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.7352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsHHH) + group(COCsFO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]OOC=[C]OF(6649)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u0 p2 c0 {1,S} {7,S}
5  C u0 p0 c0 {2,S} {7,D} {8,S}
6  C u1 p0 c0 {3,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (340.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,235.037],'cm^-1')),
        HinderedRotor(inertia=(0.0301984,'amu*angstrom^2'), symmetry=1, barrier=(51.9714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316272,'amu*angstrom^2'), symmetry=1, barrier=(12.466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.444301,'amu*angstrom^2'), symmetry=1, barrier=(17.3971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31634,'amu*angstrom^2'), symmetry=1, barrier=(51.9692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23582,0.0645356,-8.56852e-05,5.96716e-08,-1.66421e-11,41004.5,26.6009], Tmin=(100,'K'), Tmax=(874.323,'K')), NASAPolynomial(coeffs=[10.8799,0.0204138,-9.98869e-06,1.95285e-09,-1.38155e-13,39318.1,-18.6309], Tmin=(874.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CsJOOC) + radical(C=CJO)"""),
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
    E0 = (-134.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (240.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (204.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-126.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (57.2008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-46.7745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-52.7067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-92.0452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-129.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (257.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (26.8145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (76.4277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (63.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (34.3279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (527.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-126.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (181.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (22.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (371.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (213.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (38.3585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (412.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['CH2O(20)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]O[O](231)', '[CH]C(=O)F(323)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(18)', '[O]O[CH]C(=O)F(1948)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['O=C(F)C1COO1(6176)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['[CH2]OOC1O[C]1F(6639)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(192.058,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 190.6 to 192.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['F[C]1[CH]OOCO1(6640)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['[O]C1(F)[CH]OOC1(6641)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.93454e+07,'s^-1'), n=0.812267, Ea=(82.1501,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][O](2599)', 'O=CC(=O)F(335)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.01583e+10,'m^3/(mol*s)'), n=-1.23787, Ea=(159.383,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.08456738906112354, var=37.926887421795335, Tref=1000.0, N=24, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2O(20)', '[O][CH]C(=O)F(509)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(165.044,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2]OOC=C=O(6642)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(21.7215,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][O](2599)', '[O][CH]C(=O)F(509)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(9.60291,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OOC(F)[C]=O(6630)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(165.885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['O=[C][CH]OOCF(6643)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(198.406,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['[CH2]OC([O])(F)C=O(6173)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(169.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O(6)', '[CH2]OOC=[C]F(6644)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['FC1=COOCO1(6645)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Ypri_rad_out] + [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['[CH2][O](2599)', 'FC1=COO1(285)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(316.066,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OO[CH]C(=O)F(6174)'],
    products = ['[O][C](F)C1COO1(6646)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(157.099,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 156.2 to 157.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[CH2]OO[C]=C=O(6635)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(249.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OO[C]=C(O)F(6647)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['COO[C]C(=O)F(6648)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]OOC=[C]OF(6649)'],
    products = ['[CH2]OO[CH]C(=O)F(6174)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(33.7202,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1778',
    isomers = [
        '[CH2]OO[CH]C(=O)F(6174)',
    ],
    reactants = [
        ('CH2O(20)', 'O=CC(=O)F(335)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1778',
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

