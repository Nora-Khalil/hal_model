species(
    label = '[O]C(F)(C=O)O[CH]C(=O)F(1910)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u1 p0 c0 {3,S} {10,S} {12,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-673.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.506113,0.107266,-0.000180219,1.51805e-07,-4.98286e-11,-80843.8,33.5963], Tmin=(100,'K'), Tmax=(813.198,'K')), NASAPolynomial(coeffs=[14.6065,0.0248137,-1.31611e-05,2.57807e-09,-1.7915e-13,-83033.4,-34.5387], Tmin=(813.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-673.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(COCsFO) + radical(C=OCOJ) + radical(CCsJOCs)"""),
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
    label = '[O]C(F)O[CH]C(=O)F(1895)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
7  C u1 p0 c0 {3,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-548.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,3025,407.5,1350,352.5,611,648,830,1210,1753,180.12,181.036,496.441,500.468],'cm^-1')),
        HinderedRotor(inertia=(0.0521494,'amu*angstrom^2'), symmetry=1, barrier=(9.18547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26726,'amu*angstrom^2'), symmetry=1, barrier=(9.23813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136271,'amu*angstrom^2'), symmetry=1, barrier=(23.8875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3929.02,'J/mol'), sigma=(6.055,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=613.70 K, Pc=40.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00632,0.0671353,-8.42873e-05,5.2251e-08,-1.27013e-11,-65827.9,26.452], Tmin=(100,'K'), Tmax=(1007.45,'K')), NASAPolynomial(coeffs=[13.9799,0.015625,-7.59382e-06,1.50037e-09,-1.07521e-13,-68442,-36.2343], Tmin=(1007.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(CsFHOO) + group(COCsFO) + radical(O2sj(Cs-F1sO2sH)) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C(F)[CH]OC(=O)[CH]OF(2344)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u0 p2 c0 {7,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {5,D} {9,S}
8  C u1 p0 c0 {3,S} {10,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {11,S}
10 C u0 p0 c0 {1,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-462.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3000,3050,390,425,1340,1360,335,370,611,648,830,1210,1753,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0679011,0.0929464,-0.000128506,9.05014e-08,-2.54076e-11,-55478.2,32.4735], Tmin=(100,'K'), Tmax=(869.188,'K')), NASAPolynomial(coeffs=[14.4508,0.0267564,-1.42791e-05,2.88939e-09,-2.08205e-13,-57978.5,-34.8992], Tmin=(869.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-462.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(COCsFO) + radical(CCsJOC(O)) + radical(OCJC=O)"""),
)

species(
    label = 'O=C(F)[CH]OO[CH]C(=O)F(1912)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {3,S} {9,S} {11,S}
8  C u1 p0 c0 {4,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-586.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3000,3050,390,425,1340,1360,335,370,542,680,621,675,735,925,1110,1310,1739,1767,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4051.1,'J/mol'), sigma=(6.06009,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=632.77 K, Pc=41.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0462439,0.0951889,-0.000143018,1.16144e-07,-3.82908e-11,-70451.1,32.6507], Tmin=(100,'K'), Tmax=(739.265,'K')), NASAPolynomial(coeffs=[11.4016,0.0337448,-1.83401e-05,3.70471e-09,-2.65279e-13,-72129.9,-18.7016], Tmin=(739.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-586.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(COCsFO) + radical(OCJC=O) + radical(OCJC=O)"""),
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
    label = '[O]C([O])(F)C=O(1992)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-272.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.620008,'amu*angstrom^2'), symmetry=1, barrier=(14.2552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9391,0.0520955,-9.61727e-05,8.85868e-08,-3.0716e-11,-32742.1,19.1588], Tmin=(100,'K'), Tmax=(871.373,'K')), NASAPolynomial(coeffs=[6.59299,0.0158185,-8.05218e-06,1.52656e-09,-1.02948e-13,-32987,0.596067], Tmin=(871.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCOJ)"""),
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
    label = 'O=C[C](F)O[CH]C(=O)F(2345)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {1,S} {3,S} {8,S}
7  C u1 p0 c0 {3,S} {9,S} {10,S}
8  C u0 p0 c0 {4,D} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-569.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,611,648,830,1210,1753,243.643,243.68,243.691,695.249],'cm^-1')),
        HinderedRotor(inertia=(1.23607,'amu*angstrom^2'), symmetry=1, barrier=(52.0716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.405083,'amu*angstrom^2'), symmetry=1, barrier=(17.0575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284079,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23643,'amu*angstrom^2'), symmetry=1, barrier=(52.072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721521,0.0779985,-0.00011018,8.28115e-08,-2.50977e-11,-68383.6,27.8968], Tmin=(100,'K'), Tmax=(804.27,'K')), NASAPolynomial(coeffs=[11.1393,0.0261821,-1.35331e-05,2.69356e-09,-1.9185e-13,-70059.2,-20.0931], Tmin=(804.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-569.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=CC1(F)OC(C(=O)F)O1(2346)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-990.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43563,0.044361,2.06039e-05,-1.5831e-07,1.61322e-10,-119068,25.8771], Tmin=(100,'K'), Tmax=(389.565,'K')), NASAPolynomial(coeffs=[3.97122,0.0420932,-2.26425e-05,4.64966e-09,-3.39186e-13,-119290,18.6017], Tmin=(389.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-990.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(COCsFO) + ring(O2s-Cs-O2s-Cs(F))"""),
)

species(
    label = '[O]C(F)(C=O)OC1O[C]1F(2347)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
9  C u1 p0 c0 {2,S} {4,S} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-523.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.322775,0.0880666,-0.000125407,9.50716e-08,-2.91428e-11,-62821.1,30.1328], Tmin=(100,'K'), Tmax=(794.26,'K')), NASAPolynomial(coeffs=[11.8294,0.0301182,-1.59695e-05,3.2154e-09,-2.30537e-13,-64649,-22.7296], Tmin=(794.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-523.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + ring(Cs(O2)-O2s-Cs(F)) + radical(C=OCOJ) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'O=C(F)[CH]OC1(F)[CH]OO1(2348)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u1 p0 c0 {5,S} {7,S} {11,S}
9  C u1 p0 c0 {3,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-416.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.128012,0.0915605,-0.000130768,9.59533e-08,-2.80294e-11,-49952.4,30.7174], Tmin=(100,'K'), Tmax=(837.247,'K')), NASAPolynomial(coeffs=[13.8004,0.0262349,-1.37231e-05,2.74817e-09,-1.96509e-13,-52241.6,-32.8145], Tmin=(837.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-416.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCsJOO) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1(F)[CH]OC(C(=O)F)O1(2349)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-701.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.341726,0.0833204,-0.000108234,6.56637e-08,-1.44707e-11,-84246.9,28.5086], Tmin=(100,'K'), Tmax=(1324.58,'K')), NASAPolynomial(coeffs=[20.3487,0.00262677,3.76994e-06,-1.08886e-09,8.7231e-14,-88130.5,-71.0958], Tmin=(1324.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-701.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(COCsFO) + ring(1,3-Dioxolane) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CC1(F)O[CH][C](F)OO1(2253)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u1 p0 c0 {3,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {5,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-493.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07222,0.0930502,-0.000111527,6.3178e-08,-1.30708e-11,-59107.8,29.79], Tmin=(100,'K'), Tmax=(1408.67,'K')), NASAPolynomial(coeffs=[22.8747,0.00529705,2.95036e-06,-9.54863e-10,7.7762e-14,-63894.5,-86.9886], Tmin=(1408.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-OdCsH) + ring(124trioxane) + radical(CCsJOCs) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1OC1(F)O[CH]C(=O)F(2350)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u0 p0 c0 {3,S} {5,S} {7,S} {11,S}
9  C u1 p0 c0 {4,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-614.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532506,0.0834889,-0.000111161,7.86148e-08,-2.26633e-11,-73841.4,30.6251], Tmin=(100,'K'), Tmax=(839.314,'K')), NASAPolynomial(coeffs=[11.658,0.0304658,-1.63973e-05,3.34284e-09,-2.42124e-13,-75709,-21.1001], Tmin=(839.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + ring(Cs(F)(O2)-O2s-Cs) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1C(C(=O)F)OC1([O])F(2351)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-593.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.810122,0.064251,-5.68255e-05,2.4235e-08,-4.08368e-12,-71264.4,32.0258], Tmin=(100,'K'), Tmax=(1420.09,'K')), NASAPolynomial(coeffs=[16.9183,0.0188791,-8.90103e-06,1.73692e-09,-1.23034e-13,-75839.4,-51.3366], Tmin=(1420.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-593.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFOO) + group(COCsFO) + ring(O2s-Cs-Cs-Cs(F)) + radical(CC(C)OJ) + radical(O2sj(Cs-F1sO2sCs))"""),
)

species(
    label = '[O]C1(F)[CH]OC(F)(C=O)O1(2322)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {11,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-681.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0559331,0.0823907,-9.31887e-05,3.40987e-08,3.49836e-12,-81854.5,23.9505], Tmin=(100,'K'), Tmax=(797.422,'K')), NASAPolynomial(coeffs=[22.0745,0.00127949,4.20049e-06,-1.18434e-09,9.55279e-14,-86298.9,-83.1407], Tmin=(797.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-681.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(1,3-Dioxolane) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOCs)"""),
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
    label = 'O=CC(=O)O[CH]C(=O)F(2352)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {2,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,D} {8,S}
8  C u0 p0 c0 {4,D} {7,S} {11,S}
9  C u0 p0 c0 {1,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-628.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,611,648,830,1210,1753,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950329,0.073293,-9.13753e-05,5.92895e-08,-1.56885e-11,-75469.4,27.5507], Tmin=(100,'K'), Tmax=(909.791,'K')), NASAPolynomial(coeffs=[11.548,0.0266977,-1.45497e-05,2.99231e-09,-2.18239e-13,-77397.7,-22.5747], Tmin=(909.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-628.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)O2s) + group(COCsFO) + group(Cds-O2d(Cds-O2d)H) + radical(CCsJOC(O))"""),
)

species(
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C(F)[CH]OC(=O)F(2174)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {3,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-756.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,482,664,788,1296,1923,180,180,180,993.435],'cm^-1')),
        HinderedRotor(inertia=(0.437388,'amu*angstrom^2'), symmetry=1, barrier=(19.3387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841242,'amu*angstrom^2'), symmetry=1, barrier=(19.3418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0276088,'amu*angstrom^2'), symmetry=1, barrier=(19.3403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27907,0.0603605,-7.18957e-05,4.10935e-08,-9.19399e-12,-90933.4,24.5936], Tmin=(100,'K'), Tmax=(1090.49,'K')), NASAPolynomial(coeffs=[14.1284,0.0132274,-7.06182e-06,1.45698e-09,-1.07003e-13,-93735.8,-38.5102], Tmin=(1090.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(COFOO) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C(F)(C=O)OC=C=O(2353)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {9,D} {11,S}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-472.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20038,'amu*angstrom^2'), symmetry=1, barrier=(27.599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19968,'amu*angstrom^2'), symmetry=1, barrier=(27.5831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19957,'amu*angstrom^2'), symmetry=1, barrier=(27.5806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0984668,0.0880495,-0.000118173,7.47047e-08,-1.80254e-11,-56713.1,27.9172], Tmin=(100,'K'), Tmax=(1028.65,'K')), NASAPolynomial(coeffs=[20.025,0.00979827,-4.06628e-06,7.52811e-10,-5.25249e-14,-60853.2,-69.7353], Tmin=(1028.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=OCOJ)"""),
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
    label = 'O=[C]C(=O)O[CH]C(=O)F(2354)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {2,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,D} {9,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-469.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,1855,455,950,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83783,0.0743832,-9.95488e-05,6.6191e-08,-1.74634e-11,-56397.8,28.7256], Tmin=(100,'K'), Tmax=(923.999,'K')), NASAPolynomial(coeffs=[13.4477,0.0197951,-1.09321e-05,2.25418e-09,-1.64517e-13,-58728.1,-31.1132], Tmin=(923.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)O2s) + group(COCsFO) + group(Cds-O2d(Cds-O2d)H) + radical(CCsJOC(O)) + radical(OC=OCJ=O)"""),
)

species(
    label = 'O=[C]C(O)(F)O[CH]C(=O)F(2355)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {12,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u1 p0 c0 {3,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-762.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.64726,0.109951,-0.000183724,1.52109e-07,-4.91274e-11,-91550,34.0027], Tmin=(100,'K'), Tmax=(792.979,'K')), NASAPolynomial(coeffs=[15.9054,0.0229173,-1.23997e-05,2.44995e-09,-1.71238e-13,-94063.9,-41.3137], Tmin=(792.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-762.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(COCsFO) + radical(CCsJOCs) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C(F)([C]=O)OCC(=O)F(2356)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-699.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.146821,0.100562,-0.000170642,1.50162e-07,-5.17304e-11,-83959.2,33.5377], Tmin=(100,'K'), Tmax=(808.493,'K')), NASAPolynomial(coeffs=[11.4961,0.0301398,-1.62046e-05,3.20482e-09,-2.24339e-13,-85422.8,-17.5664], Tmin=(808.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(COCsFO) + radical(C=OCOJ) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=C(OF)O[CH]C(=O)F(2357)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {3,S} {10,S} {11,S}
9  C u0 p0 c0 {5,S} {7,D} {12,S}
10 C u0 p0 c0 {1,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-378.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,611,648,830,1210,1753,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03299,0.107455,-0.000147756,9.32695e-08,-2.23132e-11,-45377.6,32.8819], Tmin=(100,'K'), Tmax=(1040.11,'K')), NASAPolynomial(coeffs=[25.0273,0.00723301,-3.21878e-06,6.26738e-10,-4.54389e-14,-50798.7,-93.8687], Tmin=(1040.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-378.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(COCsFO) + radical(C=COJ) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C=C([O])OC(F)C(=O)F(2358)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {10,D}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-817.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.706198,0.0932147,-0.000106517,5.3332e-08,-9.29336e-12,-98173.2,31.2996], Tmin=(100,'K'), Tmax=(1023.07,'K')), NASAPolynomial(coeffs=[25.4608,0.00657323,-2.45502e-06,4.89056e-10,-3.79572e-14,-104347,-99.5448], Tmin=(1023.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-817.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(COCsFO) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=COC(F)(C=O)OF(2359)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-382.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00186,0.102669,-0.000132316,7.81959e-08,-1.74745e-11,-45790.8,34.2834], Tmin=(100,'K'), Tmax=(1114.49,'K')), NASAPolynomial(coeffs=[25.9085,0.00608598,-2.32422e-06,4.37807e-10,-3.20401e-14,-51789.1,-98.4608], Tmin=(1114.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-382.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)(C=O)OC(F)[C]=O(2329)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-670.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,355,410,600,1181,1341,1420,3056,2782.5,750,1395,475,1775,1000,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.64416,'amu*angstrom^2'), symmetry=1, barrier=(14.8105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645085,'amu*angstrom^2'), symmetry=1, barrier=(14.8318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64458,'amu*angstrom^2'), symmetry=1, barrier=(37.8121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644485,'amu*angstrom^2'), symmetry=1, barrier=(14.818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.530121,0.112693,-0.000206855,1.87374e-07,-6.43998e-11,-80435.2,33.11], Tmin=(100,'K'), Tmax=(857.683,'K')), NASAPolynomial(coeffs=[11.5322,0.0300032,-1.60082e-05,3.09551e-09,-2.11465e-13,-81532.1,-17.5639], Tmin=(857.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-670.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(COj(Cs-F1sO2sH)(O2d))"""),
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
    label = 'O=C(F)[CH]O(1552)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-466.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,634.702],'cm^-1')),
        HinderedRotor(inertia=(0.178809,'amu*angstrom^2'), symmetry=1, barrier=(51.1156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178808,'amu*angstrom^2'), symmetry=1, barrier=(51.1155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9486,0.00327696,8.62766e-05,-1.78856e-07,1.15095e-10,-56066.9,9.87278], Tmin=(10,'K'), Tmax=(486.336,'K')), NASAPolynomial(coeffs=[2.12805,0.0278262,-1.89743e-05,5.9068e-09,-6.93439e-13,-56003.1,16.1793], Tmin=(486.336,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-466.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)[CH]O""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)(C=O)C([O])(F)C=O(1909)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {7,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-613.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31031,'amu*angstrom^2'), symmetry=1, barrier=(30.1267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31031,'amu*angstrom^2'), symmetry=1, barrier=(30.1265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30865,'amu*angstrom^2'), symmetry=1, barrier=(30.0885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.202884,0.10088,-0.00016716,1.43629e-07,-4.85395e-11,-73687.1,31.7586], Tmin=(100,'K'), Tmax=(800.937,'K')), NASAPolynomial(coeffs=[12.3798,0.0289184,-1.53083e-05,3.01553e-09,-2.10914e-13,-75410.1,-24.3263], Tmin=(800.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(F)(C=O)OC=[C]F(2360)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
7  C u0 p0 c0 {5,D} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {9,D} {11,S}
9  C u1 p0 c0 {2,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-351.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14908,'amu*angstrom^2'), symmetry=1, barrier=(26.4196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15005,'amu*angstrom^2'), symmetry=1, barrier=(26.4419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15011,'amu*angstrom^2'), symmetry=1, barrier=(26.4433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0160438,0.0932422,-0.000141348,1.05158e-07,-3.04547e-11,-42135.6,28.7903], Tmin=(100,'K'), Tmax=(851.804,'K')), NASAPolynomial(coeffs=[15.9499,0.0182685,-9.32299e-06,1.82918e-09,-1.28683e-13,-44855.6,-45.6752], Tmin=(851.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-351.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=OCOJ) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'O=CC1(F)OC=C(F)OO1(2259)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u0 p0 c0 {3,S} {9,D} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,D}
10 C u0 p0 c0 {6,D} {7,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-684.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275832,0.0718391,-7.01958e-05,3.35168e-08,-6.21562e-12,-82142.1,25.6071], Tmin=(100,'K'), Tmax=(1321.07,'K')), NASAPolynomial(coeffs=[18.9893,0.015178,-5.8609e-06,1.05112e-09,-7.18528e-14,-87086.5,-69.8855], Tmin=(1321.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCFOO) + group(Cds-CdsOsH) + group(CdCFO) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
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
    label = '[O][C](F)C1OC(F)(C=O)O1(2361)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {5,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-619.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881633,0.0795063,-0.000134253,1.32025e-07,-5.09423e-11,-74440,31.9321], Tmin=(100,'K'), Tmax=(799.543,'K')), NASAPolynomial(coeffs=[3.99961,0.0419318,-2.25321e-05,4.49511e-09,-3.17469e-13,-74236.1,21.9798], Tmin=(799.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFOO) + group(CsCFHO) + group(Cds-OdCsH) + ring(O2s-Cs-O2s-Cs(F)) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1OC(F)=COC1([O])F(2362)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {4,S} {6,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u0 p0 c0 {2,S} {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-580.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515345,0.074069,-8.38182e-05,4.95193e-08,-1.15206e-11,-69693.9,21.3062], Tmin=(100,'K'), Tmax=(1054.44,'K')), NASAPolynomial(coeffs=[14.6525,0.0204399,-7.52776e-06,1.28488e-09,-8.45496e-14,-72675.3,-47.6469], Tmin=(1054.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(CdCFO) + ring(23dihydro14dioxin) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCOJ) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C(F)C=O(506)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-324.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([241,322,589,815,1096,1220,1302,2892,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.586986,'amu*angstrom^2'), symmetry=1, barrier=(13.496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82948,0.0149886,6.15997e-05,-2.22156e-07,2.15915e-10,-39081.1,10.2785], Tmin=(10,'K'), Tmax=(359.038,'K')), NASAPolynomial(coeffs=[4.10807,0.0226824,-1.56543e-05,5.05177e-09,-6.15171e-13,-39170.7,8.25062], Tmin=(359.038,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-324.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]C(F)CDO""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)(C=O)O[C]=C=O(2338)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u1 p0 c0 {2,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-233.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03523,'amu*angstrom^2'), symmetry=1, barrier=(23.8021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03432,'amu*angstrom^2'), symmetry=1, barrier=(23.7811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03418,'amu*angstrom^2'), symmetry=1, barrier=(23.7778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.137605,0.0894761,-0.000143835,1.11922e-07,-3.36364e-11,-27893.2,30.2003], Tmin=(100,'K'), Tmax=(823.788,'K')), NASAPolynomial(coeffs=[15.5956,0.0144309,-7.21209e-06,1.37634e-09,-9.42827e-14,-30440.5,-41.382], Tmin=(823.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)(C=O)O[C]=C(O)F(2363)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {9,S} {12,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-572.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,293,496,537,1218,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.541782,0.106105,-0.000165853,1.27122e-07,-3.79313e-11,-68749.3,33.8572], Tmin=(100,'K'), Tmax=(827.018,'K')), NASAPolynomial(coeffs=[17.1398,0.0205942,-1.07733e-05,2.12429e-09,-1.49392e-13,-71674.2,-48.09], Tmin=(827.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)=[C]OC(O)(F)C=O(2364)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {7,S} {12,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {6,S} {10,D}
10 C u1 p0 c0 {3,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-675.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.226759,0.10027,-0.000151996,1.16063e-07,-3.50132e-11,-81061.8,33.3622], Tmin=(100,'K'), Tmax=(813.453,'K')), NASAPolynomial(coeffs=[15.0176,0.0253098,-1.37724e-05,2.78319e-09,-1.99321e-13,-83542,-37.0356], Tmin=(813.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-675.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)([C]=O)OC=C(O)F(2365)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {9,S} {12,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {9,D} {11,S}
9  C u0 p0 c0 {2,S} {4,S} {8,D}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,326,540,652,719,1357,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16817,'amu*angstrom^2'), symmetry=1, barrier=(26.8586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16757,'amu*angstrom^2'), symmetry=1, barrier=(26.8448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16703,'amu*angstrom^2'), symmetry=1, barrier=(26.8323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16908,'amu*angstrom^2'), symmetry=1, barrier=(26.8794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.66621,0.104314,-0.000147139,9.79999e-08,-2.50816e-11,-78972.4,31.8528], Tmin=(100,'K'), Tmax=(966.635,'K')), NASAPolynomial(coeffs=[21.1451,0.0140592,-7.08671e-06,1.41053e-09,-1.01245e-13,-83189.2,-72.6341], Tmin=(966.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=OCOJ) + radical(CsCJ=O)"""),
)

species(
    label = 'O=CC(=O)O[CH][C](F)OF(2366)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {2,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u1 p0 c0 {3,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {5,D} {10,S}
9  C u1 p0 c0 {1,S} {4,S} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-364.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,395,473,707,1436,2782.5,750,1395,475,1775,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.580993,0.112424,-0.000199611,1.80075e-07,-6.30652e-11,-43655.6,33.7345], Tmin=(100,'K'), Tmax=(811.516,'K')), NASAPolynomial(coeffs=[11.9757,0.032439,-1.83244e-05,3.67291e-09,-2.58334e-13,-45097.8,-20.5514], Tmin=(811.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(CCsJOC(O)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)(C=O)OC=[C]OF(2367)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {2,S} {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-256.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.630072,0.103748,-0.000149738,1.02899e-07,-2.71618e-11,-30709.7,34.6218], Tmin=(100,'K'), Tmax=(938.391,'K')), NASAPolynomial(coeffs=[20.3622,0.0142666,-6.70269e-06,1.28136e-09,-8.95136e-14,-34649.5,-65.3186], Tmin=(938.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=CJO)"""),
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
    E0 = (-234.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (77.5269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (190.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (41.9024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (161.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (112.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-226.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-51.1199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (22.4431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-114.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-54.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-132.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-154.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-130.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-183.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-88.0616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-161.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (78.9298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (14.8693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-6.84426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-159.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-128.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (142.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-73.1815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (140.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-44.8041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-170.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-41.8566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (330.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-227.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (122.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-112.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-112.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-56.0425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (215.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (16.3014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-88.9721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-134.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (152.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (234.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=CC(=O)F(335)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)O[CH]C(=O)F(1895)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.81927e-16,'m^3/(mol*s)'), n=6.80628, Ea=(305.597,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C(F)[CH]OC(=O)[CH]OF(2344)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(213.916,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C(F)[CH]OO[CH]C(=O)F(1912)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.25418e+20,'s^-1'), n=-1.9758, Ea=(189.903,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(=O)F(323)', '[O]C([O])(F)C=O(1992)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', 'O=C[C](F)O[CH]C(=O)F(2345)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=CC1(F)OC(C(=O)F)O1(2346)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C(F)(C=O)OC1O[C]1F(2347)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=C(F)[CH]OC1(F)[CH]OO1(2348)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(257.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 255.4 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C1(F)[CH]OC(C(=O)F)O1(2349)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.46737e+16,'s^-1'), n=-1.32849, Ea=(120.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHCO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=CC1(F)O[CH][C](F)OO1(2253)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.32511e+10,'s^-1'), n=0.274726, Ea=(180.356,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 177.8 to 180.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C1OC1(F)O[CH]C(=O)F(2350)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C1C(C(=O)F)OC1([O])F(2351)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(79.9378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Exocyclic
Ea raised from 74.2 to 79.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C1(F)[CH]OC(F)(C=O)O1(2322)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.16887e+10,'s^-1'), n=0.492453, Ea=(104.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_O] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CC(=O)F(335)', '[O][CH]C(=O)F(509)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(75.5886,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=CC(=O)O[CH]C(=O)F(2352)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(28.5222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HCO(15)', 'O=C(F)[CH]OC(=O)F(2174)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(124.435,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F(37)', '[O]C(F)(C=O)OC=C=O(2353)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.9385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][CH]C(=O)F(509)', '[O][CH]C(=O)F(509)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(4.94065,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'O=[C]C(=O)O[CH]C(=O)F(2354)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(305.218,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=[C]C(O)(F)O[CH]C(=O)F(2355)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C(F)([C]=O)OCC(=O)F(2356)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.78e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/OneDe;XH_out] for rate rule [R4H_SSS_OCs;C_rad_out_H/OneDe;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C=C(OF)O[CH]C(=O)F(2357)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(82.7094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C=C([O])OC(F)C(=O)F(2358)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(161.394,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=COC(F)(C=O)OF(2359)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(83.4726,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(F)(C=O)OC(F)[C]=O(2329)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(186.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=[C]C(=O)F(1543)', 'O=C(F)[CH]O(1552)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.52781e-06,'m^3/(mol*s)'), n=3.38614, Ea=(181.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OR] for rate rule [Cdd_Cd;H_OCdpri]
Euclidian distance = 2.23606797749979
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(F)(C=O)C([O])(F)C=O(1909)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.87873e+12,'s^-1'), n=0.324012, Ea=(133.121,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O(6)', '[O]C(F)(C=O)OC=[C]F(2360)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=CC1(F)OC=C(F)OO1(2259)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C1(F)[CH]OOC(F)=CO1(2286)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(357.448,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 349.7 to 357.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O][C](F)C1OC(F)(C=O)O1(2361)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C1OC(F)=COC1([O])F(2362)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_O] for rate rule [R7_SMSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.8284271247461903
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['O=[C]C(=O)F(1543)', '[O]C(F)C=O(506)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.33333e+07,'s^-1'), n=1.2, Ea=(178.533,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_N-5CO->O"""),
)

reaction(
    label = 'reaction35',
    reactants = ['HF(38)', '[O]C(F)(C=O)O[C]=C=O(2338)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(290.474,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C(F)(C=O)O[C]=C(O)F(2363)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    products = ['[O]C(F)=[C]OC(O)(F)C=O(2364)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C(F)([C]=O)OC=C(O)F(2365)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(16140,'s^-1'), n=1.92259, Ea=(84.7802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;CO_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O=CC(=O)O[CH][C](F)OF(2366)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(78.0414,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]C(F)(C=O)OC=[C]OF(2367)'],
    products = ['[O]C(F)(C=O)O[CH]C(=O)F(1910)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(52.0458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #355',
    isomers = [
        '[O]C(F)(C=O)O[CH]C(=O)F(1910)',
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
    label = 'PDepNetwork #355',
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

