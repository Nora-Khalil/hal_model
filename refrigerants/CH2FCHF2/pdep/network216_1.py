species(
    label = 'F[C](F)CC1[C](F)C2OOC12F(1014)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {7,S} {12,S} {15,S} {16,S}
11 C u1 p0 c0 {2,S} {7,S} {9,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-442.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,2750,2850,1437.5,1250,1305,750,350,212,367,445,1450,190,488,555,1236,1407,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.212558,0.0946781,-9.55049e-05,5.18623e-08,-1.21422e-11,-53122.3,28.8909], Tmin=(100,'K'), Tmax=(985.92,'K')), NASAPolynomial(coeffs=[11.2985,0.049701,-2.70755e-05,5.59098e-09,-4.09123e-13,-55308.3,-24.435], Tmin=(985.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-442.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'FC1=CC2(F)OOC12(939)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-248.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,323,467,575,827,1418,556.295,556.354,556.401,556.506,556.667,556.693,556.734,556.988,556.993],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3642.96,'J/mol'), sigma=(5.89185,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.02 K, Pc=40.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37505,0.0476133,-2.48125e-05,-2.36513e-09,3.49961e-12,-29797.1,16.5191], Tmin=(100,'K'), Tmax=(1189.18,'K')), NASAPolynomial(coeffs=[15.9615,0.016946,-9.33452e-06,1.96667e-09,-1.45909e-13,-34567,-61.8485], Tmin=(1189.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'F[C](F)CC1[C](F)C(F)C12OO2(1131)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {7,S} {12,S} {15,S} {16,S}
11 C u1 p0 c0 {2,S} {7,S} {9,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-435.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,212,367,445,1450,190,488,555,1236,1407,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12894,0.10689,-0.000112334,5.72342e-08,-1.14102e-11,-52168.9,35.2484], Tmin=(100,'K'), Tmax=(1221.26,'K')), NASAPolynomial(coeffs=[23.8845,0.0249628,-1.17069e-05,2.30263e-09,-1.65208e-13,-58278.4,-90.4264], Tmin=(1221.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s1_3_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)C[CH]C1(F)C2OOC21F(1132)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {11,S} {12,S} {14,S} {15,S}
11 C u1 p0 c0 {7,S} {10,S} {16,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-437.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0842281,0.101369,-0.000115673,7.24571e-08,-1.93353e-11,-52494.4,31.8346], Tmin=(100,'K'), Tmax=(882.91,'K')), NASAPolynomial(coeffs=[11.4147,0.0492734,-2.71666e-05,5.62791e-09,-4.12377e-13,-54524.9,-22.2093], Tmin=(882.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-437.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_3_4_ane) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)CC1(F)[CH]C2(F)OOC12(1013)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {7,S} {12,S} {14,S} {15,S}
11 C u1 p0 c0 {7,S} {9,S} {16,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-448.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0223281,0.100306,-0.000109539,6.5089e-08,-1.65562e-11,-53853.1,29.1069], Tmin=(100,'K'), Tmax=(919.545,'K')), NASAPolynomial(coeffs=[11.5124,0.0501287,-2.76855e-05,5.74367e-09,-4.21305e-13,-55974.3,-25.5738], Tmin=(919.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-448.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)CC1C2(F)[CH]OOC12F(1133)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {12,S} {14,S} {15,S}
11 C u1 p0 c0 {6,S} {8,S} {16,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-495.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.122683,0.0948321,-9.42878e-05,4.97134e-08,-1.11299e-11,-59462.9,29.7066], Tmin=(100,'K'), Tmax=(1036.93,'K')), NASAPolynomial(coeffs=[12.5624,0.0468457,-2.48721e-05,5.08471e-09,-3.70103e-13,-62042.8,-30.759], Tmin=(1036.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-495.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + polycyclic(s2_3_5_ane) + radical(CCsJOOC) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-465.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.141016,0.100934,-0.000104417,5.5609e-08,-1.23631e-11,-55799.9,28.3539], Tmin=(100,'K'), Tmax=(1053.25,'K')), NASAPolynomial(coeffs=[14.6145,0.0448961,-2.46116e-05,5.09555e-09,-3.73346e-13,-58908.2,-43.5986], Tmin=(1053.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-465.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'F[C]F(156)',
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
    label = '[CH2]C1[C](F)C2OOC12F(1134)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  C u1 p0 c0 {6,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-4.07728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01003,0.0633488,-4.71787e-05,1.65746e-08,-2.33571e-12,-381.432,22.4476], Tmin=(100,'K'), Tmax=(1630.68,'K')), NASAPolynomial(coeffs=[16.3084,0.0258221,-1.26592e-05,2.46197e-09,-1.72091e-13,-5370.76,-58.8391], Tmin=(1630.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.07728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Isobutyl)"""),
)

species(
    label = 'FC12OOC1C1(F)C2CC1(F)F(1135)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {10,S}
10 C u0 p0 c0 {6,S} {8,S} {9,S} {14,S}
11 C u0 p0 c0 {7,S} {12,S} {15,S} {16,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-690.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381215,0.0859338,-6.80088e-05,2.52059e-08,-3.81356e-12,-82871.3,19.4095], Tmin=(100,'K'), Tmax=(1488.8,'K')), NASAPolynomial(coeffs=[17.578,0.0397309,-2.14585e-05,4.36127e-09,-3.13331e-13,-87991.8,-70.3988], Tmin=(1488.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-690.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) - ring(Cyclobutane)"""),
)

species(
    label = 'FC1=C(CC(F)F)C2(F)OOC12(1136)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {10,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {16,S}
11 C u0 p0 c0 {7,S} {9,S} {12,D}
12 C u0 p0 c0 {4,S} {8,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-720.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.31667,0.0957232,-8.72705e-05,3.813e-08,-6.69671e-12,-86519.1,27.968], Tmin=(100,'K'), Tmax=(1338.75,'K')), NASAPolynomial(coeffs=[20.3308,0.0340319,-1.81491e-05,3.70942e-09,-2.69019e-13,-92047.6,-77.668], Tmin=(1338.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'FC(F)=CC1C(F)C2OOC12F(1137)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {10,S} {14,S}
10 C u0 p0 c0 {2,S} {7,S} {9,S} {15,S}
11 C u0 p0 c0 {7,S} {12,D} {16,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-687.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.557063,0.0893901,-7.21485e-05,2.69808e-08,-3.95748e-12,-82500.5,28.2267], Tmin=(100,'K'), Tmax=(1617.09,'K')), NASAPolynomial(coeffs=[25.2776,0.0254855,-1.28706e-05,2.54253e-09,-1.79322e-13,-90855.8,-108.827], Tmin=(1617.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-687.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFF) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'FC1=C2OOC2(F)C1CC(F)F(1015)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {16,S}
11 C u0 p0 c0 {6,S} {8,S} {12,D}
12 C u0 p0 c0 {4,S} {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-650.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.274459,0.0934064,-9.90546e-05,5.91888e-08,-1.54778e-11,-78123.6,26.2138], Tmin=(100,'K'), Tmax=(886.823,'K')), NASAPolynomial(coeffs=[9.58584,0.0514077,-2.80172e-05,5.78688e-09,-4.23594e-13,-79775.1,-17.5899], Tmin=(886.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane)"""),
)

species(
    label = 'F[C](F)CC1[C](F)OOC=C1F(1138)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {7,S} {12,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u1 p0 c0 {2,S} {5,S} {7,S}
11 C u0 p0 c0 {6,S} {9,D} {16,S}
12 C u1 p0 c0 {3,S} {4,S} {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-504.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,395,473,707,1436,190,488,555,1236,1407,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.441118,0.104152,-0.000132623,9.36501e-08,-2.70037e-11,-60561.6,32.9406], Tmin=(100,'K'), Tmax=(840.914,'K')), NASAPolynomial(coeffs=[12.9331,0.0405347,-1.91468e-05,3.68817e-09,-2.58697e-13,-62811,-29.2652], Tmin=(840.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-504.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(CsCsF1sO2s) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C](F)CC=C(F)C1OO[C]1F(1139)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {13,S}
8  C u0 p0 c0 {10,S} {12,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 C u1 p0 c0 {2,S} {6,S} {7,S}
12 C u1 p0 c0 {3,S} {4,S} {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-452.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,3010,987.5,1337.5,450,1655,395,473,707,1436,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.142001,0.0966052,-0.000106559,6.13805e-08,-1.44345e-11,-54332.4,39.7779], Tmin=(100,'K'), Tmax=(1017.61,'K')), NASAPolynomial(coeffs=[15.0483,0.0368949,-1.85426e-05,3.718e-09,-2.68219e-13,-57424,-33.7717], Tmin=(1017.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHO) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
)

species(
    label = '[O]OC1(F)C=C(F)C1C[C](F)F(1140)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {7,S} {12,S} {14,S} {15,S}
10 C u0 p0 c0 {2,S} {7,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {16,S}
12 C u1 p0 c0 {3,S} {4,S} {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-531.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,190,488,555,1236,1407,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.527895,0.106848,-0.000131797,8.65061e-08,-2.31039e-11,-63774.1,34.6063], Tmin=(100,'K'), Tmax=(904.484,'K')), NASAPolynomial(coeffs=[14.7443,0.0393062,-1.97834e-05,3.94239e-09,-2.82702e-13,-66536.7,-37.5397], Tmin=(904.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-531.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2][C](F)F(163)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25818,0.0186259,-1.56413e-05,8.22305e-09,-2.04579e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.66,'K')), NASAPolynomial(coeffs=[4.47009,0.0131647,-6.4128e-06,1.29206e-09,-9.37476e-14,-13305.9,7.85289], Tmin=(887.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
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
    label = 'F[C](F)CC1=C(F)C2OOC12F(1141)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {6,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {10,S} {12,S} {14,S} {15,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {2,S} {8,S} {10,D}
12 C u1 p0 c0 {3,S} {4,S} {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-519.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (183.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.225436,0.0926469,-8.56348e-05,3.74539e-08,-6.53916e-12,-62314.2,29.5752], Tmin=(100,'K'), Tmax=(1352.9,'K')), NASAPolynomial(coeffs=[20.8501,0.0303349,-1.65478e-05,3.41003e-09,-2.48279e-13,-68016.8,-78.4722], Tmin=(1352.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-519.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)CC1C(F)=C2OOC21F(1142)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {12,S} {14,S} {15,S}
10 C u0 p0 c0 {6,S} {8,S} {11,D}
11 C u0 p0 c0 {2,S} {7,S} {10,D}
12 C u1 p0 c0 {3,S} {4,S} {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
"""),
    E0 = (-450.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (183.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404065,0.0932531,-9.70939e-05,2.16723e-08,2.93724e-11,-54013.1,27.211], Tmin=(100,'K'), Tmax=(500.431,'K')), NASAPolynomial(coeffs=[7.98821,0.0515222,-2.86301e-05,5.89574e-09,-4.2854e-13,-55008.7,-6.49098], Tmin=(500.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(CdCsCdF) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C]1[CH]C2(F)OOC12(943)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {5,S} {8,S} {10,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (24.0771,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,3150,900,1100,212,367,445,1450,792.942,792.942,792.942,792.942,792.943,792.943,792.943,792.943,792.943,792.943],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01799,0.0472587,-3.27286e-05,9.66921e-09,-1.11055e-12,2962.68,18.0674], Tmin=(100,'K'), Tmax=(1964.1,'K')), NASAPolynomial(coeffs=[16.3603,0.0180501,-1.04219e-05,2.09781e-09,-1.46839e-13,-2671.3,-60.8074], Tmin=(1964.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.0771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'F[C]1C(C=C(F)F)C2(F)OOC12(1143)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {10,S} {14,S}
10 C u1 p0 c0 {2,S} {7,S} {9,S}
11 C u0 p0 c0 {7,S} {12,D} {15,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-496.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,212,367,445,1450,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (183.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279887,0.0845846,-7.16213e-05,2.89666e-08,-4.75903e-12,-59633.6,27.8923], Tmin=(100,'K'), Tmax=(1403.76,'K')), NASAPolynomial(coeffs=[17.7051,0.0349319,-1.85648e-05,3.7695e-09,-2.71627e-13,-64525.8,-62.0844], Tmin=(1403.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-496.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFF) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s)"""),
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
    label = 'F[C](F)CC1=C(F)C2OO[C]12(1144)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {4,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {8,S} {11,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 C u0 p0 c0 {1,S} {6,S} {8,D}
11 C u1 p0 c0 {2,S} {3,S} {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-88.1026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,190,488,555,1236,1407,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.322708,0.0770164,-6.3566e-05,2.43492e-08,-3.69679e-12,-10461.3,31.5932], Tmin=(100,'K'), Tmax=(1544.18,'K')), NASAPolynomial(coeffs=[20.2398,0.0254243,-1.34506e-05,2.71316e-09,-1.93992e-13,-16612.5,-73.1492], Tmin=(1544.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.1026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOO) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)CC1[C](F)C2=C1OO2(1145)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {6,S} {11,S} {13,S} {14,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u1 p0 c0 {1,S} {6,S} {10,S}
10 C u0 p0 c0 {5,S} {8,D} {9,S}
11 C u1 p0 c0 {2,S} {3,S} {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (21.7808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,346,659,817,1284,190,488,555,1236,1407,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582994,0.0717369,-5.53303e-05,1.96675e-08,-2.77264e-12,2744.72,31.017], Tmin=(100,'K'), Tmax=(1645.46,'K')), NASAPolynomial(coeffs=[19.5367,0.0256623,-1.33294e-05,2.65087e-09,-1.87285e-13,-3492.89,-69.8635], Tmin=(1645.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.7808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + polycyclic(s2_4_4_ene_m) + radical(CsCdCsF1s) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[C](F)C[C]1C(F)C2OOC12F(1146)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {14,S}
10 C u0 p0 c0 {11,S} {12,S} {15,S} {16,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-442.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444306,0.0815598,-6.11528e-05,2.10775e-08,-2.92398e-12,-53073,30.8453], Tmin=(100,'K'), Tmax=(1626.73,'K')), NASAPolynomial(coeffs=[19.0261,0.0358689,-1.90215e-05,3.81134e-09,-2.7048e-13,-59118.5,-67.8423], Tmin=(1626.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-442.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)CC1C(F)[C]2OOC21F(1147)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {14,S}
10 C u0 p0 c0 {7,S} {12,S} {15,S} {16,S}
11 C u1 p0 c0 {6,S} {8,S} {9,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-446.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,316,385,515,654,689,1295,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.252476,0.0943047,-9.57912e-05,5.37359e-08,-1.31909e-11,-53587.4,26.849], Tmin=(100,'K'), Tmax=(937.436,'K')), NASAPolynomial(coeffs=[10.1065,0.0522579,-2.8512e-05,5.88982e-09,-4.31125e-13,-55434.9,-20.0546], Tmin=(937.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-446.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1C([CH]C(F)F)C2(F)OOC12(1089)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {12,S} {15,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 C u1 p0 c0 {7,S} {10,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-448.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448741,0.0897612,-8.21338e-05,3.88638e-08,-7.92419e-12,-53852.7,29.223], Tmin=(100,'K'), Tmax=(1108.74,'K')), NASAPolynomial(coeffs=[11.8156,0.0487529,-2.66539e-05,5.50458e-09,-4.02283e-13,-56373.3,-26.7889], Tmin=(1108.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-448.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Cs_S)"""),
)

species(
    label = 'F[C](F)[CH]C1C(F)C2OOC12F(1148)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {10,S} {14,S}
10 C u0 p0 c0 {2,S} {7,S} {9,S} {15,S}
11 C u1 p0 c0 {7,S} {12,S} {16,S}
12 C u1 p0 c0 {3,S} {4,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-438.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,250,417,511,1155,1315,1456,3119,3025,407.5,1350,352.5,190,488,555,1236,1407,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440037,0.08884,-7.81254e-05,3.44296e-08,-6.43516e-12,-52640.9,27.8987], Tmin=(100,'K'), Tmax=(1207.17,'K')), NASAPolynomial(coeffs=[13.3055,0.0462102,-2.5155e-05,5.17665e-09,-3.77033e-13,-55747.1,-36.5922], Tmin=(1207.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-438.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1[C](CC(F)F)C2(F)OOC12(1149)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {10,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {16,S}
11 C u1 p0 c0 {7,S} {9,S} {12,S}
12 C u1 p0 c0 {4,S} {8,S} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-452.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.609092,0.0809206,-6.0841e-05,2.12013e-08,-3.00778e-12,-54292.3,31.5888], Tmin=(100,'K'), Tmax=(1571.43,'K')), NASAPolynomial(coeffs=[17.144,0.0388318,-2.06655e-05,4.15724e-09,-2.96228e-13,-59489,-55.656], Tmin=(1571.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'F[C]1[C]2OOC2(F)C1CC(F)F(1150)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {16,S}
11 C u1 p0 c0 {6,S} {8,S} {12,S}
12 C u1 p0 c0 {4,S} {7,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-456.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.136122,0.0967901,-0.00010557,6.59429e-08,-1.80897e-11,-54794.2,28.614], Tmin=(100,'K'), Tmax=(848.218,'K')), NASAPolynomial(coeffs=[9.2702,0.0537162,-2.93977e-05,6.07514e-09,-4.44693e-13,-56343.7,-13.9491], Tmin=(848.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = 'F[C](F)CC1[C]2OOC2C1(F)F(1151)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {12,S} {15,S} {16,S}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 C u1 p0 c0 {3,S} {4,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-434.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,222,329,445,522,589,1214,1475,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125343,0.0963002,-0.000106292,6.72578e-08,-1.85541e-11,-52152.7,29.3617], Tmin=(100,'K'), Tmax=(847.661,'K')), NASAPolynomial(coeffs=[9.48631,0.0521264,-2.81221e-05,5.77789e-09,-4.21579e-13,-53739.6,-14.2523], Tmin=(847.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-434.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1C(CC(F)(F)F)[C]2OOC12(1152)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {11,S} {12,S} {16,S}
10 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
11 C u1 p0 c0 {6,S} {7,S} {9,S}
12 C u1 p0 c0 {4,S} {7,S} {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
"""),
    E0 = (-459.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50584,0.0875506,-8.40775e-05,4.51321e-08,-1.07185e-11,-55180.1,29.4909], Tmin=(100,'K'), Tmax=(959.668,'K')), NASAPolynomial(coeffs=[9.25577,0.0510804,-2.70736e-05,5.53286e-09,-4.02756e-13,-56859.5,-12.3624], Tmin=(959.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFF) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
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
    E0 = (-139.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (166.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (19.9539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (19.9539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (19.9539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (105.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (332.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-131.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-76.5811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-76.5811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-115.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-71.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-25.0254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-139.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-30.8943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (10.0294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (64.5129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-28.1294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (24.9954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (217.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (161.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (225.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (96.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (34.0734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-19.9273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (24.1159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-3.90149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-17.6243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-100.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (91.7713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (93.0176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['CH2CF2(57)', 'FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)CC1[C](F)C(F)C12OO2(1131)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C](F)C[CH]C1(F)C2OOC21F(1132)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C](F)CC1(F)[CH]C2(F)OOC12(1013)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C](F)CC1C2(F)[CH]OOC12F(1133)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]F(156)', '[CH2]C1[C](F)C2OOC12F(1134)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['FC12OOC1C1(F)C2CC1(F)F(1135)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['FC1=C(CC(F)F)C2(F)OOC12(1136)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['FC(F)=CC1C(F)C2OOC12F(1137)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['FC1=C2OOC2(F)C1CC(F)F(1015)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C](F)CC1[C](F)OOC=C1F(1138)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c6_gamma;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C](F)CC=C(F)C1OO[C]1F(1139)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC1(F)C=C(F)C1C[C](F)F(1140)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.10864e+11,'s^-1'), n=0.310877, Ea=(88.8126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 88.0 to 88.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C](F)F(163)', 'FC1=CC2(F)OOC12(939)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(23.9976,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'F[C](F)CC1=C(F)C2OOC12F(1141)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(205.9,'m^3/(mol*s)'), n=1.596, Ea=(14.8388,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'F[C](F)CC1C(F)=C2OOC21F(1142)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.296163,'m^3/(mol*s)'), n=2.44476, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_Ext-9R!H-R_Ext-8R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_Ext-9R!H-R_Ext-8R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2CF2(57)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(6.65128,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'F[C]1C(C=C(F)F)C2(F)OOC12(1143)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(301,'m^3/(mol*s)'), n=1.6, Ea=(7.34249,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS_Ext-5R!H-R_Ext-6R!H-R_Ext-7R!H-R_N-1COS-inRing',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS_Ext-5R!H-R_Ext-6R!H-R_Ext-7R!H-R_N-1COS-inRing"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C](F)F(163)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'F[C](F)CC1=C(F)C2OO[C]12(1144)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(228.218,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'F[C](F)CC1[C](F)C2=C1OO2(1145)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(181.713,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CF2(43)', '[CH2]C1[C](F)C2OOC12F(1134)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(2.01137,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C](F)C[C]1C(F)C2OOC12F(1146)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.93363e+09,'s^-1'), n=1.033, Ea=(174.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_Cs2] for rate rule [R2H_S_cy4;C_rad_out_noH;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C](F)CC1C(F)[C]2OOC21F(1147)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.26888e+09,'s^-1'), n=0.972714, Ea=(123.905,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_NDMustO;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C]1C([CH]C(F)F)C2(F)OOC12(1089)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.23544e+08,'s^-1'), n=1.37575, Ea=(164.097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] + [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[C](F)[CH]C1C(F)C2OOC12F(1148)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(17992.9,'s^-1'), n=2.5106, Ea=(132.012,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_23cy4;Y_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C]1[C](CC(F)F)C2(F)OOC12(1149)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.50974e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C]1[C]2OOC2(F)C1CC(F)F(1150)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(10179.1,'s^-1'), n=1.84089, Ea=(39.9141,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_single;XH_out] for rate rule [R5H_CCC;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[C](F)CC1[C]2OOC2C1(F)F(1151)'],
    products = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(223.716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C](F)CC1[C](F)C2OOC12F(1014)'],
    products = ['F[C]1C(CC(F)(F)F)[C]2OOC12(1152)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(232.999,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #216',
    isomers = [
        'F[C](F)CC1[C](F)C2OOC12F(1014)',
    ],
    reactants = [
        ('CH2CF2(57)', 'FC1=CC2(F)OOC12(939)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #216',
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

