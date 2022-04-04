species(
    label = 'O=C[C](F)OO[CH]C(=O)F(1913)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {3,S} {9,S} {11,S}
8  C u1 p0 c0 {1,S} {4,S} {10,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-560.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,280,501,1494,1531,611,648,830,1210,1753,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.12626,'amu*angstrom^2'), symmetry=1, barrier=(48.8869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12599,'amu*angstrom^2'), symmetry=1, barrier=(48.8807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12609,'amu*angstrom^2'), symmetry=1, barrier=(48.883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12693,'amu*angstrom^2'), symmetry=1, barrier=(48.9024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429865,'amu*angstrom^2'), symmetry=1, barrier=(48.8835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767327,0.0826676,-9.75387e-05,3.70121e-08,1.48088e-11,-67343.5,29.8178], Tmin=(100,'K'), Tmax=(527.422,'K')), NASAPolynomial(coeffs=[9.10531,0.0374561,-2.02178e-05,4.07296e-09,-2.9115e-13,-68473.8,-7.45054], Tmin=(527.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(COCsFO) + radical(CsCOF1sO2s) + radical(OCJC=O)"""),
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
    label = '[O]O[C](F)C=O(1934)',
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
    label = 'O=C[C]F(284)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u2 p0 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (22.7053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,163,1167],'cm^-1')),
        HinderedRotor(inertia=(0.0337629,'amu*angstrom^2'), symmetry=1, barrier=(13.6228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34132,0.0165481,-1.72264e-05,1.26789e-08,-4.5452e-12,2752.64,10.5202], Tmin=(100,'K'), Tmax=(638.16,'K')), NASAPolynomial(coeffs=[4.07864,0.0119262,-6.36171e-06,1.32811e-09,-9.81565e-14,2658.54,7.2943], Tmin=(638.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.7053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
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
    label = 'O=CC1(F)OOC1C(=O)F(1916)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-744.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689393,0.0629903,-5.27952e-05,2.07242e-08,-3.17509e-12,-89439.7,31.648], Tmin=(100,'K'), Tmax=(1567.13,'K')), NASAPolynomial(coeffs=[19.2164,0.015701,-7.53134e-06,1.46862e-09,-1.03286e-13,-95246.6,-66.0574], Tmin=(1567.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-744.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(COCsFO) + ring(Cs-Cs(F)-O2s-O2s)"""),
)

species(
    label = 'O=C=C(F)OOCC(=O)F(1918)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-615.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82037,0.0817158,-9.8852e-05,3.549e-08,1.82426e-11,-73970.1,30.9165], Tmin=(100,'K'), Tmax=(522.544,'K')), NASAPolynomial(coeffs=[9.4829,0.0346246,-1.88434e-05,3.80086e-09,-2.71376e-13,-75137.8,-7.76357], Tmin=(522.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C[C](F)OOC1O[C]1F(2271)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {1,S} {3,S} {7,S}
9  C u1 p0 c0 {2,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-368.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00646,0.0730717,-7.72727e-05,4.39071e-08,-1.05657e-11,-44247.5,28.0578], Tmin=(100,'K'), Tmax=(974.114,'K')), NASAPolynomial(coeffs=[10.2696,0.0350346,-1.87008e-05,3.82145e-09,-2.7803e-13,-46052.2,-16.3885], Tmin=(974.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C(F)[CH]OOC1(F)[CH]O1(2272)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {11,S}
9  C u1 p0 c0 {5,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-414.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519705,0.0827289,-9.80276e-05,6.00749e-08,-1.50092e-11,-49760.1,29.7571], Tmin=(100,'K'), Tmax=(960.893,'K')), NASAPolynomial(coeffs=[13.1088,0.0303221,-1.62165e-05,3.31329e-09,-2.41013e-13,-52179.4,-30.4755], Tmin=(960.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + ring(Cs(F)(O2)-O2s-Cs) + radical(Csj(Cs-F1sO2sO2s)(O2s-Cs)(H)_ring) + radical(OCJC=O)"""),
)

species(
    label = 'O=CC1(F)OO[CH][C](F)O1(2273)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-487.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.560731,0.085421,-9.75584e-05,5.45715e-08,-1.13364e-11,-58409.6,28.8081], Tmin=(100,'K'), Tmax=(1377.41,'K')), NASAPolynomial(coeffs=[20.1948,0.0104932,-3.97845e-09,-3.68541e-10,3.71426e-14,-62737.3,-72.9256], Tmin=(1377.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-487.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(124trioxane) + radical(CsCsF1sO2s) + radical(CCsJOOC)"""),
)

species(
    label = 'O=C(F)C1O[CH][C](F)OO1(2274)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u1 p0 c0 {3,S} {9,S} {12,S}
9  C u1 p0 c0 {1,S} {5,S} {8,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-513.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.977157,0.0874881,-0.000100308,5.46976e-08,-1.08514e-11,-61521,32.6228], Tmin=(100,'K'), Tmax=(1487.76,'K')), NASAPolynomial(coeffs=[21.6599,0.0058547,2.93902e-06,-9.51781e-10,7.67001e-14,-65957.9,-77.8552], Tmin=(1487.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-513.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(COCsFO) + ring(124trioxane) + radical(CCsJOCs) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)[CH]OOC1(F)C=O(2275)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-456.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.293276,0.0755063,-7.46056e-05,3.1869e-08,-3.22038e-12,-54774.8,27.5759], Tmin=(100,'K'), Tmax=(917.44,'K')), NASAPolynomial(coeffs=[17.5553,0.0166123,-5.07605e-06,7.90954e-10,-5.07492e-14,-58631,-57.9701], Tmin=(917.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C1[C](F)OOC1C(=O)F(2276)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {12,S}
9  C u1 p0 c0 {1,S} {4,S} {8,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-468.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.509765,0.0758919,-7.58101e-05,3.64338e-08,-6.49347e-12,-56142.5,34.2194], Tmin=(100,'K'), Tmax=(1617,'K')), NASAPolynomial(coeffs=[20.6489,0.00911588,-4.74642e-07,-1.46885e-10,1.57613e-14,-61097.9,-72.1915], Tmin=(1617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(COCsFO) + ring(12dioxolane) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C=COO[C](F)C=O(2277)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {1,S} {2,S} {8,S}
7  C u0 p0 c0 {3,S} {9,D} {10,S}
8  C u0 p0 c0 {4,D} {6,S} {11,S}
9  C u0 p0 c0 {5,D} {7,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-263.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,280,501,1494,1531,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(1.91433,'amu*angstrom^2'), symmetry=1, barrier=(44.0142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91435,'amu*angstrom^2'), symmetry=1, barrier=(44.0146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871565,'amu*angstrom^2'), symmetry=1, barrier=(20.039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91448,'amu*angstrom^2'), symmetry=1, barrier=(44.0176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.860341,0.0749649,-0.000100991,7.36262e-08,-2.18826e-11,-31536,27.8468], Tmin=(100,'K'), Tmax=(815.967,'K')), NASAPolynomial(coeffs=[10.47,0.0278574,-1.4394e-05,2.87434e-09,-2.05515e-13,-33104.2,-16.5598], Tmin=(815.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCOF1sO2s)"""),
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
    label = 'O=C=C(F)OO[CH]C(=O)F(2278)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {3,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-477.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,611,648,830,1210,1753,197,221,431,657,2120,512.5,787.5,180,180,1413.32],'cm^-1')),
        HinderedRotor(inertia=(0.391437,'amu*angstrom^2'), symmetry=1, barrier=(8.9999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0331813,'amu*angstrom^2'), symmetry=1, barrier=(47.067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04581,'amu*angstrom^2'), symmetry=1, barrier=(47.0371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0331522,'amu*angstrom^2'), symmetry=1, barrier=(47.0455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.185789,0.0940919,-0.000162446,1.48909e-07,-5.35622e-11,-57289.6,32.8325], Tmin=(100,'K'), Tmax=(798.925,'K')), NASAPolynomial(coeffs=[9.07756,0.0335235,-1.85939e-05,3.72614e-09,-2.63041e-13,-58198.2,-4.86342], Tmin=(798.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-477.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(OCJC=O)"""),
)

species(
    label = 'O=C[C]F-2(1596)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (35.6539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,262,1290],'cm^-1')),
        HinderedRotor(inertia=(0.407026,'amu*angstrom^2'), symmetry=1, barrier=(9.35834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79853,0.0176142,-2.7066e-05,3.1192e-08,-1.61007e-11,4288.74,9.0865], Tmin=(10,'K'), Tmax=(518.444,'K')), NASAPolynomial(coeffs=[4.38817,0.0119962,-7.71993e-06,2.33921e-09,-2.7033e-13,4241.97,6.76768], Tmin=(518.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(35.6539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]C(F)OO[CH]C(=O)F(2279)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
8  C u1 p0 c0 {4,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-541.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,355,410,600,1181,1341,1420,3056,3025,407.5,1350,352.5,611,648,830,1210,1753,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.977147,'amu*angstrom^2'), symmetry=1, barrier=(22.4665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11704,'amu*angstrom^2'), symmetry=1, barrier=(48.6748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11681,'amu*angstrom^2'), symmetry=1, barrier=(48.6697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335334,'amu*angstrom^2'), symmetry=1, barrier=(7.70998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.117,'amu*angstrom^2'), symmetry=1, barrier=(48.6739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.413581,0.106695,-0.000180427,1.57711e-07,-5.40511e-11,-64975.7,33.4056], Tmin=(100,'K'), Tmax=(804.995,'K')), NASAPolynomial(coeffs=[12.3164,0.0311797,-1.68698e-05,3.33983e-09,-2.33931e-13,-66628,-22.7804], Tmin=(804.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-541.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-OdCsH) + radical(OCJC=O) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C][C](F)OOCC(=O)F(2280)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u1 p0 c0 {2,S} {4,S} {10,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-540.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,486,617,768,1157,1926,280,501,1494,1531,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.30021,'amu*angstrom^2'), symmetry=1, barrier=(52.8864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30715,'amu*angstrom^2'), symmetry=1, barrier=(53.0459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30739,'amu*angstrom^2'), symmetry=1, barrier=(53.0515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217878,'amu*angstrom^2'), symmetry=1, barrier=(5.00944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784625,'amu*angstrom^2'), symmetry=1, barrier=(18.0401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.073775,0.0993117,-0.000166674,1.48293e-07,-5.18133e-11,-64899,32.6458], Tmin=(100,'K'), Tmax=(809.569,'K')), NASAPolynomial(coeffs=[10.2918,0.0335389,-1.78346e-05,3.51404e-09,-2.45685e-13,-66100.2,-12.2259], Tmin=(809.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-540.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C]C(F)OO[C](F)C=O(2256)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
8  C u1 p0 c0 {2,S} {4,S} {9,S}
9  C u0 p0 c0 {5,D} {8,S} {12,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-515.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,355,410,600,1181,1341,1420,3056,280,501,1494,1531,2782.5,750,1395,475,1775,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.91507,'amu*angstrom^2'), symmetry=1, barrier=(44.0312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92736,'amu*angstrom^2'), symmetry=1, barrier=(44.3138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93108,'amu*angstrom^2'), symmetry=1, barrier=(44.3993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.557682,'amu*angstrom^2'), symmetry=1, barrier=(12.8222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91897,'amu*angstrom^2'), symmetry=1, barrier=(44.1209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.129346,0.101151,-0.000171496,1.54202e-07,-5.4364e-11,-61849.5,32.0431], Tmin=(100,'K'), Tmax=(809.546,'K')), NASAPolynomial(coeffs=[10.0019,0.0349021,-1.87447e-05,3.70567e-09,-2.59487e-13,-62959.3,-11.4172], Tmin=(809.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C][CH]OOC(F)(F)C=O(2281)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {11,S}
9  C u1 p0 c0 {4,S} {10,S} {12,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-514.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,251,367,519,700,855,1175,1303,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.91176,'amu*angstrom^2'), symmetry=1, barrier=(43.9551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91409,'amu*angstrom^2'), symmetry=1, barrier=(44.0087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91175,'amu*angstrom^2'), symmetry=1, barrier=(43.9548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155626,'amu*angstrom^2'), symmetry=1, barrier=(3.57814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47329,'amu*angstrom^2'), symmetry=1, barrier=(33.8738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0372384,0.0971919,-0.000150714,1.25328e-07,-4.21204e-11,-61704.7,32.7503], Tmin=(100,'K'), Tmax=(726.783,'K')), NASAPolynomial(coeffs=[11.6331,0.0329501,-1.81014e-05,3.66291e-09,-2.62119e-13,-63400.7,-19.8258], Tmin=(726.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-514.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C(O[CH]C(=O)F)C(=O)F(1907)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
8  C u1 p0 c0 {3,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-704.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0610863,0.093451,-0.000143074,1.12695e-07,-3.52087e-11,-84656.2,35.4887], Tmin=(100,'K'), Tmax=(785.424,'K')), NASAPolynomial(coeffs=[13.3856,0.0255897,-1.34673e-05,2.68106e-09,-1.90105e-13,-86749.2,-25.5755], Tmin=(785.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(COCsFO) + radical(C=OCOJ) + radical(CCsJOCs)"""),
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
    label = '[CH]=C(F)OO[CH]C(=O)F(2282)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {7,D}
6  C u1 p0 c0 {3,S} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {5,D} {6,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u1 p0 c0 {8,D} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-184.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,611,648,830,1210,1753,293,496,537,1218,3120,650,792.5,1650,180,1731.62],'cm^-1')),
        HinderedRotor(inertia=(1.93634,'amu*angstrom^2'), symmetry=1, barrier=(44.5203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93622,'amu*angstrom^2'), symmetry=1, barrier=(44.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266887,'amu*angstrom^2'), symmetry=1, barrier=(6.13626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93684,'amu*angstrom^2'), symmetry=1, barrier=(44.5318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.322043,0.0899866,-0.000151557,1.36872e-07,-4.89095e-11,-22081.4,30.9576], Tmin=(100,'K'), Tmax=(788.001,'K')), NASAPolynomial(coeffs=[9.27996,0.0320735,-1.7633e-05,3.53353e-09,-2.49901e-13,-23106.9,-7.67386], Tmin=(788.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + group(COCsFO) + group(Cds-CdsHH) + radical(OCJC=O) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'O=C(F)C1OC=C(F)OO1(2283)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {9,D} {12,S}
9  C u0 p0 c0 {1,S} {5,S} {8,D}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-704.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728013,0.0622415,-4.56912e-05,9.00091e-09,2.32384e-12,-84571.1,27.1471], Tmin=(100,'K'), Tmax=(1023.89,'K')), NASAPolynomial(coeffs=[16.9933,0.0171761,-6.74001e-06,1.26454e-09,-9.06468e-14,-88870.4,-56.4371], Tmin=(1023.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(CdCFO) + group(COCsFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
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
    label = 'F[C]1[CH]OOC(F)=COO1(2284)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u1 p0 c0 {3,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  C u0 p0 c0 {1,S} {4,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-50.3096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08201,0.0525578,-3.08236e-05,7.23373e-09,-6.28273e-13,-5992.16,23.4222], Tmin=(100,'K'), Tmax=(2772.47,'K')), NASAPolynomial(coeffs=[35.6993,0.00407079,-4.59841e-06,9.29527e-10,-5.998e-14,-24638.4,-173.053], Tmin=(2772.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.3096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFHO) + group(CdCFO) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(CCsJOOC) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O][CH]C1(F)OOC1C(=O)F(2285)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-417.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83155,0.0768117,-9.73433e-05,6.75089e-08,-1.94128e-11,-50108,31.3426], Tmin=(100,'K'), Tmax=(835.369,'K')), NASAPolynomial(coeffs=[10.1618,0.0321352,-1.71208e-05,3.48645e-09,-2.52669e-13,-51666.8,-11.9918], Tmin=(835.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(COCsFO) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C1(F)[CH]OOC(F)=CO1(2286)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
8  C u1 p0 c0 {4,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u0 p0 c0 {2,S} {5,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-316.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247719,0.039482,9.84509e-05,-1.99242e-07,9.14639e-11,-37831.8,27.8857], Tmin=(100,'K'), Tmax=(910.694,'K')), NASAPolynomial(coeffs=[40.6022,-0.0245355,1.73942e-05,-3.37926e-09,2.18161e-13,-49877.3,-188.805], Tmin=(910.694,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(CdCFO) + ring(Cycloheptane) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    label = 'O=C=[C]OO[CH]C(=O)F(2287)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {2,S} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {4,D} {6,S}
8  C u1 p0 c0 {3,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-49.4417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,611,648,830,1210,1753,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.7981,'amu*angstrom^2'), symmetry=1, barrier=(41.3419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79595,'amu*angstrom^2'), symmetry=1, barrier=(41.2924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489168,'amu*angstrom^2'), symmetry=1, barrier=(11.2469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0389102,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524756,0.0856621,-0.000150287,1.35977e-07,-4.76352e-11,-5830.22,32.5019], Tmin=(100,'K'), Tmax=(825.87,'K')), NASAPolynomial(coeffs=[9.29187,0.0270961,-1.46672e-05,2.89035e-09,-2.01289e-13,-6729.14,-4.7922], Tmin=(825.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.4417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(OCJC=O) + radical(C=CJO)"""),
)

species(
    label = 'O=C(F)[CH]OOC(F)=[C]O(2288)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {10,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {3,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-389.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3025,407.5,1350,352.5,293,496,537,1218,611,648,830,1210,1753,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.236797,0.103373,-0.000176684,1.57872e-07,-5.53276e-11,-46682.8,36.4466], Tmin=(100,'K'), Tmax=(802.34,'K')), NASAPolynomial(coeffs=[10.9741,0.0330087,-1.80769e-05,3.6013e-09,-2.53253e-13,-48015.9,-12.2672], Tmin=(802.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COCsFO) + group(Cds-CdsOsH) + radical(OCJC=O) + radical(C=CJO)"""),
)

species(
    label = 'O=C(F)[CH]OO[C]=COF(2289)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {2,S} {9,S}
6  O u0 p2 c0 {8,D}
7  C u1 p0 c0 {3,S} {8,S} {11,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-73.1114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,3025,407.5,1350,352.5,611,648,830,1210,1753,3010,987.5,1337.5,450,1655,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1701,0.0989422,-0.000152061,1.20607e-07,-3.79878e-11,-8649.8,36.6706], Tmin=(100,'K'), Tmax=(778.935,'K')), NASAPolynomial(coeffs=[13.7395,0.0275136,-1.45106e-05,2.88179e-09,-2.03886e-13,-10816.7,-26.96], Tmin=(778.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(COCsFO) + group(Cds-CdsOsH) + radical(OCJC=O) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OOC(F)C(=O)F(2290)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-456.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,3010,987.5,1337.5,450,1655,1685,370,202.324,203.211,205.847],'cm^-1')),
        HinderedRotor(inertia=(0.827411,'amu*angstrom^2'), symmetry=1, barrier=(24.5683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822272,'amu*angstrom^2'), symmetry=1, barrier=(24.4984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19755,'amu*angstrom^2'), symmetry=1, barrier=(35.5065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84879,'amu*angstrom^2'), symmetry=1, barrier=(24.5854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.148529,0.0870977,-0.000106854,6.43775e-08,-1.52828e-11,-54778.2,34.8383], Tmin=(100,'K'), Tmax=(1027.2,'K')), NASAPolynomial(coeffs=[16.9811,0.0215505,-1.11371e-05,2.25611e-09,-1.63808e-13,-58236.3,-46.8212], Tmin=(1027.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(COCsFO) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]OOC(F)=COF(2291)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {11,S}
9  C u1 p0 c0 {4,S} {10,S} {12,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-62.6038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,326,540,652,719,1357,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1855,455,950,180,2418.14],'cm^-1')),
        HinderedRotor(inertia=(0.0786899,'amu*angstrom^2'), symmetry=1, barrier=(1.80924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89257,'amu*angstrom^2'), symmetry=1, barrier=(43.5139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88569,'amu*angstrom^2'), symmetry=1, barrier=(43.3557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88664,'amu*angstrom^2'), symmetry=1, barrier=(43.3776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88583,'amu*angstrom^2'), symmetry=1, barrier=(43.359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.16248,0.103676,-0.000182377,1.69982e-07,-6.17827e-11,-7391.37,35.2693], Tmin=(100,'K'), Tmax=(805.636,'K')), NASAPolynomial(coeffs=[8.78907,0.038194,-2.12872e-05,4.26707e-09,-3.00958e-13,-8150.99,-1.74499], Tmin=(805.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.6038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C(F)(C=O)O[C](F)C=O(1911)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u1 p0 c0 {2,S} {3,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-689.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0627585,0.097645,-0.000168627,1.54588e-07,-5.50738e-11,-82777.9,31.5212], Tmin=(100,'K'), Tmax=(823.11,'K')), NASAPolynomial(coeffs=[8.71625,0.0355371,-1.88967e-05,3.71443e-09,-2.58846e-13,-83523,-4.41538], Tmin=(823.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-689.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C[C](F)OOC=[C]F(2292)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u1 p0 c0 {1,S} {3,S} {8,S}
7  C u0 p0 c0 {4,S} {9,D} {10,S}
8  C u0 p0 c0 {5,D} {6,S} {11,S}
9  C u1 p0 c0 {2,S} {7,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-141.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,280,501,1494,1531,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(2.19663,'amu*angstrom^2'), symmetry=1, barrier=(50.5049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333297,'amu*angstrom^2'), symmetry=1, barrier=(7.66315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19573,'amu*angstrom^2'), symmetry=1, barrier=(50.4842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1958,'amu*angstrom^2'), symmetry=1, barrier=(50.4857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.569758,0.0850028,-0.000143356,1.32518e-07,-4.82108e-11,-16942.6,30.0297], Tmin=(100,'K'), Tmax=(803.287,'K')), NASAPolynomial(coeffs=[7.45497,0.0344335,-1.85191e-05,3.67644e-09,-2.58482e-13,-17523.4,1.59093], Tmin=(803.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsCOF1sO2s) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'O=CC1(F)OOC=C(F)O1(2293)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-684.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275832,0.0718391,-7.01958e-05,3.35168e-08,-6.21562e-12,-82142.1,25.6071], Tmin=(100,'K'), Tmax=(1321.07,'K')), NASAPolynomial(coeffs=[18.9893,0.015178,-5.8609e-06,1.05112e-09,-7.18528e-14,-87086.5,-69.8855], Tmin=(1321.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFOO) + group(CdCFO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
)

species(
    label = 'O=C=C(F)OOC=C(O)F(2294)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-519.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215315,0.0932282,-0.000137186,9.75478e-08,-2.32033e-11,-62364.5,32.444], Tmin=(100,'K'), Tmax=(605.951,'K')), NASAPolynomial(coeffs=[11.8916,0.0303103,-1.64871e-05,3.31742e-09,-2.3633e-13,-64039.5,-20.1829], Tmin=(605.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-519.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = '[O]C1OC(F)=COO[C]1F(2295)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
8  C u1 p0 c0 {1,S} {4,S} {7,S}
9  C u0 p0 c0 {2,S} {3,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-332.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433154,0.0387376,8.88167e-05,-1.80366e-07,8.22083e-11,-39842.1,30.1855], Tmin=(100,'K'), Tmax=(917.045,'K')), NASAPolynomial(coeffs=[37.446,-0.0190746,1.38696e-05,-2.65197e-09,1.66862e-13,-50988.1,-168.934], Tmin=(917.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CdCFO) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCOJ) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'O=C=[C]OO[C](F)C=O(2264)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {1,S} {2,S} {7,S}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u1 p0 c0 {3,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-23.3572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,280,501,1494,1531,2782.5,750,1395,475,1775,1000,1685,370,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.15568,'amu*angstrom^2'), symmetry=1, barrier=(3.5794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.37706,'amu*angstrom^2'), symmetry=1, barrier=(54.6532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.38025,'amu*angstrom^2'), symmetry=1, barrier=(54.7267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.65056,'amu*angstrom^2'), symmetry=1, barrier=(106.925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813228,0.0800637,-0.000141138,1.32143e-07,-4.77897e-11,-2704.11,31.1245], Tmin=(100,'K'), Tmax=(827.016,'K')), NASAPolynomial(coeffs=[6.96712,0.0308371,-1.65533e-05,3.25892e-09,-2.27077e-13,-3056.42,6.62821], Tmin=(827.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.3572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCOF1sO2s) + radical(C=CJO)"""),
)

species(
    label = 'O=C[C](F)OO[C]=C(O)F(2296)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u0 p0 c0 {2,S} {5,S} {10,D}
9  C u0 p0 c0 {6,D} {7,S} {11,S}
10 C u1 p0 c0 {4,S} {8,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-363.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,280,501,1494,1531,293,496,537,1218,2782.5,750,1395,475,1775,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.87699,'amu*angstrom^2'), symmetry=1, barrier=(43.1558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00141618,'amu*angstrom^2'), symmetry=1, barrier=(6.63902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87744,'amu*angstrom^2'), symmetry=1, barrier=(43.166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87713,'amu*angstrom^2'), symmetry=1, barrier=(43.1588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87643,'amu*angstrom^2'), symmetry=1, barrier=(43.1429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0473642,0.0978306,-0.000167756,1.54363e-07,-5.56394e-11,-43556.5,35.0843], Tmin=(100,'K'), Tmax=(806.812,'K')), NASAPolynomial(coeffs=[8.66034,0.0367298,-1.99509e-05,3.96694e-09,-2.78792e-13,-44347.5,-0.908186], Tmin=(806.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(C=CJO)"""),
)

species(
    label = 'O=CC(F)OO[C]C(=O)F(2297)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
8  C u0 p0 c0 {5,D} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {6,D} {10,S}
10 C u2 p0 c0 {4,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-377.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000,486,617,768,1157,1926,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101483,0.0934975,-0.000126529,8.64624e-08,-2.36979e-11,-45271.6,30.8841], Tmin=(100,'K'), Tmax=(886.039,'K')), NASAPolynomial(coeffs=[14.5559,0.0282429,-1.60573e-05,3.34142e-09,-2.4479e-13,-47833,-37.1012], Tmin=(886.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(COCsFO) + radical(CH2_triplet)"""),
)

species(
    label = 'O=[C][C](F)OOC=C(O)F(2298)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u1 p0 c0 {2,S} {4,S} {10,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-444.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,280,501,1494,1531,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.92985,'amu*angstrom^2'), symmetry=1, barrier=(44.3711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652692,'amu*angstrom^2'), symmetry=1, barrier=(15.0067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52588,'amu*angstrom^2'), symmetry=1, barrier=(35.083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652256,'amu*angstrom^2'), symmetry=1, barrier=(14.9966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93079,'amu*angstrom^2'), symmetry=1, barrier=(44.3927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.411305,0.10636,-0.000180334,1.56206e-07,-5.27136e-11,-53304.8,33.2817], Tmin=(100,'K'), Tmax=(821.722,'K')), NASAPolynomial(coeffs=[12.7729,0.0290941,-1.54001e-05,3.01165e-09,-2.09037e-13,-55029.7,-25.047], Tmin=(821.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-444.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=C[C](F)OOC=[C]OF(2299)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {2,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {1,S} {3,S} {9,S}
8  C u0 p0 c0 {4,S} {10,D} {11,S}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-47.0268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,280,501,1494,1531,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1685,370,250.558,251.63],'cm^-1')),
        HinderedRotor(inertia=(0.096441,'amu*angstrom^2'), symmetry=1, barrier=(4.3172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20769,'amu*angstrom^2'), symmetry=1, barrier=(53.8438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19694,'amu*angstrom^2'), symmetry=1, barrier=(53.8385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533282,'amu*angstrom^2'), symmetry=1, barrier=(23.7039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20306,'amu*angstrom^2'), symmetry=1, barrier=(53.8437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0520568,0.0942333,-0.000146601,1.22567e-07,-4.11665e-11,-5520.88,35.5246], Tmin=(100,'K'), Tmax=(752.985,'K')), NASAPolynomial(coeffs=[11.5219,0.0310575,-1.62766e-05,3.22085e-09,-2.27153e-13,-7184.53,-16.1338], Tmin=(752.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.0268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CsCOF1sO2s) + radical(C=CJO)"""),
)

species(
    label = 'FC1=COOC(F)=COO1(2300)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {3,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u0 p0 c0 {2,S} {5,S} {7,D}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-205.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16579,0.0424779,-1.35314e-05,-3.46134e-09,1.61939e-12,-24701.4,23.9842], Tmin=(100,'K'), Tmax=(1679.16,'K')), NASAPolynomial(coeffs=[15.5463,0.0271654,-1.46474e-05,2.85555e-09,-1.95601e-13,-31529.9,-54.4564], Tmin=(1679.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclooctane)"""),
)

species(
    label = 'O=[C]C(=O)F(1543)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u1 p0 c0 {3,D} {4,S}
"""),
    E0 = (-325.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([233,496,705,1150,2014,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.06462,'amu*angstrom^2'), symmetry=1, barrier=(24.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0185,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90371,0.00648742,6.66758e-05,-1.78696e-07,1.33913e-10,-39126.2,10.0186], Tmin=(10,'K'), Tmax=(477.098,'K')), NASAPolynomial(coeffs=[4.98548,0.0142778,-1.08248e-05,3.66779e-09,-4.58866e-13,-39421.3,3.58924], Tmin=(477.098,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-325.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OD[C]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]CC(=O)F(1554)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-363.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86108,0.0121233,3.70793e-05,-8.49933e-08,5.10894e-11,-43750.9,11.7734], Tmin=(10,'K'), Tmax=(581.186,'K')), NASAPolynomial(coeffs=[3.98792,0.0215968,-1.40743e-05,4.31471e-09,-5.0294e-13,-43940.4,9.72707], Tmin=(581.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-363.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]CC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)[C]OOC(F)=CO(2301)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {11,S}
9  C u0 p0 c0 {2,S} {6,D} {10,S}
10 C u2 p0 c0 {4,S} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-306.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,326,540,652,719,1357,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.263382,0.0998089,-0.000139867,9.53472e-08,-2.55181e-11,-36723.4,32.4123], Tmin=(100,'K'), Tmax=(915.568,'K')), NASAPolynomial(coeffs=[17.5309,0.022064,-1.24884e-05,2.59258e-09,-1.89675e-13,-39981.6,-51.8642], Tmin=(915.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(COCsFO) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C=[C]OOC=C(F)OF(2302)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (26.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,2995,3025,975,1000,1300,1375,400,500,1630,1680,326,540,652,719,1357,1685,370,284.288,284.288,284.293],'cm^-1')),
        HinderedRotor(inertia=(0.217347,'amu*angstrom^2'), symmetry=1, barrier=(12.4653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698982,'amu*angstrom^2'), symmetry=1, barrier=(4.00877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570976,'amu*angstrom^2'), symmetry=1, barrier=(32.7466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08582,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0203867,0.0948852,-0.000141637,1.09553e-07,-3.38264e-11,3283.04,37.0684], Tmin=(100,'K'), Tmax=(792.376,'K')), NASAPolynomial(coeffs=[13.2722,0.0279916,-1.50104e-05,3.02038e-09,-2.15929e-13,1182.87,-23.7804], Tmin=(792.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
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
    E0 = (-243.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (122.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (124.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-235.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-204.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-51.6777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-60.2801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-170.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-168.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-139.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-151.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-235.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (148.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (51.4758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-89.3527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (137.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-73.2233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-162.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-32.4507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-21.8523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-78.7862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (375.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-236.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (53.1697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (266.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-100.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1.05415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (246.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (90.5133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (299.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (11.1418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (290.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-72.2164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (418.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-236.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-201.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-15.5779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (262.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (104.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-16.1492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-48.1484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (303.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (111.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (-157.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (35.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (366.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=CC(=O)F(335)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=O)F(323)', '[O]O[C](F)C=O(1934)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C[C]F(284)', '[O]O[CH]C(=O)F(1948)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=CC1(F)OOC1C(=O)F(1916)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=C=C(F)OOCC(=O)F(1918)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=C[C](F)OOC1O[C]1F(2271)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(192.058,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 190.6 to 192.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=C(F)[CH]OOC1(F)[CH]O1(2272)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=CC1(F)OO[CH][C](F)O1(2273)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(73.682,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 72.2 to 73.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=C(F)C1O[CH][C](F)OO1(2274)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.00901e+09,'s^-1'), n=0.463766, Ea=(75.7069,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHDe] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_csHCO]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O]C1(F)[CH]OOC1(F)C=O(2275)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.54917e+06,'s^-1'), n=1.13913, Ea=(104.219,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 102.2 to 104.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O]C1[C](F)OOC1C(=O)F(2276)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.63572e+08,'s^-1'), n=0.688067, Ea=(92.4873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_csHDe] + [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 89.0 to 92.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=CC(=O)F(335)', '[O][CH]C(=O)F(509)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(145.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O=C=COO[C](F)C=O(2277)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(21.7215,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'O=C=C(F)OO[CH]C(=O)F(2278)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C(=O)F(509)', '[O][CH]C(=O)F(509)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(22.5377,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C[C]F-2(1596)', '[O]O[CH]C(=O)F(1948)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]C(F)OO[CH]C(=O)F(2279)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.81542e+08,'s^-1'), n=1.41854, Ea=(151.197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out] for rate rule [R2H_S;CO_rad_out;Cs_H_out_OOH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C][C](F)OOCC(=O)F(2280)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(60863,'s^-1'), n=1.86284, Ea=(61.1759,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;CO_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.7320508075688772
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C]C(F)OO[C](F)C=O(2256)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(165.885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=[C][CH]OOC(F)(F)C=O(2281)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(175.269,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O]C(O[CH]C(=O)F)C(=O)F(1907)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(164.949,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(6)', '[CH]=C(F)OO[CH]C(=O)F(2282)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=C(F)C1OC=C(F)OO1(2283)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['FC1=COO1(285)', '[O][CH]C(=O)F(509)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.22709e+11,'s^-1'), n=0, Ea=(296.905,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['F[C]1[CH]OOC(F)=COO1(2284)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(510.489,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra;radadd_intra] for rate rule [R8_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic
Ea raised from 506.6 to 510.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O][CH]C1(F)OOC1C(=O)F(2285)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(143.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 141.7 to 143.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O]C1(F)[CH]OOC(F)=CO1(2286)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(244.789,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_O] + [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 238.9 to 244.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['HF(38)', 'O=C=[C]OO[CH]C(=O)F(2287)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(260.128,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=C(F)[CH]OOC(F)=[C]O(2288)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=C(F)[CH]OO[C]=COF(2289)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(55.6221,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C=[C]OOC(F)C(=O)F(2290)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(150.667,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=[C][CH]OOC(F)=COF(2291)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(35.8127,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O]C(F)(C=O)O[C](F)C=O(1911)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(171.519,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(6)', 'O=C[C](F)OOC=[C]F(2292)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=CC1(F)OOC=C(F)O1(2293)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=C=C(F)OOC=C(O)F(2294)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(41.8646,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['[O]C1OC(F)=COO[C]1F(2295)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(228.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_O] + [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 222.0 to 228.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['HF(38)', 'O=C=[C]OO[C](F)C=O(2264)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(249.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O=C[C](F)OO[C]=C(O)F(2296)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O=CC(F)OO[C]C(=O)F(2297)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O=[C][C](F)OOC=C(O)F(2298)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(364667,'s^-1'), n=1.22214, Ea=(79.2357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;Y_rad_out;XH_out] for rate rule [R7HJ_1;CO_rad_out;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O=C[C](F)OOC=[C]OF(2299)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(33.7202,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction43',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['FC1=COOC(F)=COO1(2300)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.35773e+13,'s^-1'), n=0.0154583, Ea=(354.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Ypri_rad_out] for rate rule [R8;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 348.6 to 354.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction44',
    reactants = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    products = ['O=[C]C(=O)F(1543)', '[O]CC(=O)F(1554)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(86.2167,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O=C(F)[C]OOC(F)=CO(2301)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(7.546e+09,'s^-1'), n=0.732, Ea=(25.1375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Cd_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O]C=[C]OOC=C(F)OF(2302)'],
    products = ['O=C[C](F)OO[CH]C(=O)F(1913)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(22.8626,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #358',
    isomers = [
        'O=C[C](F)OO[CH]C(=O)F(1913)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'O=CC(=O)F(335)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #358',
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

