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
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,215,315,519,588,595,1205,1248,212,367,445,1450,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.141016,0.100934,-0.000104417,5.5609e-08,-1.23631e-11,-55799.9,28.3539], Tmin=(100,'K'), Tmax=(1053.25,'K')), NASAPolynomial(coeffs=[14.6145,0.0448961,-2.46116e-05,5.09555e-09,-3.73346e-13,-58908.2,-43.5986], Tmin=(1053.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-465.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C(F)(F)C1[C](F)C(F)C12OO2(1196)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-457.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,212,367,445,1450,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5948,0.114325,-0.000124764,6.48324e-08,-1.30299e-11,-54841.3,35.1246], Tmin=(100,'K'), Tmax=(1220.76,'K')), NASAPolynomial(coeffs=[27.0444,0.0204847,-9.45811e-06,1.86276e-09,-1.34343e-13,-61833.6,-108.756], Tmin=(1220.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-457.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s1_3_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)(F)[CH]C1(F)C2OOC21F(1197)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
11 C u1 p0 c0 {7,S} {10,S} {14,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-459.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368151,0.106778,-0.000121537,7.21455e-08,-1.77808e-11,-55174.9,31.05], Tmin=(100,'K'), Tmax=(962.214,'K')), NASAPolynomial(coeffs=[14.5478,0.044771,-2.48729e-05,5.1719e-09,-3.79818e-13,-58045.3,-40.3362], Tmin=(962.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-459.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_3_4_ane) + radical(Cs_S) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {7,S} {9,S} {14,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-451.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2750,3150,900,1100,316,385,515,654,689,1295,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.317526,0.103789,-0.000108654,5.76538e-08,-1.26213e-11,-54139,28.8166], Tmin=(100,'K'), Tmax=(1076.19,'K')), NASAPolynomial(coeffs=[16.1353,0.042637,-2.34198e-05,4.8542e-09,-3.55906e-13,-57680.3,-51.7667], Tmin=(1076.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-451.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1C2(F)[CH]OOC12F(1198)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {6,S} {8,S} {14,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-517.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.272519,0.101539,-0.000104602,5.50662e-08,-1.19608e-11,-62138.6,29.3217], Tmin=(100,'K'), Tmax=(1085.2,'K')), NASAPolynomial(coeffs=[15.8515,0.0421074,-2.24557e-05,4.60214e-09,-3.35499e-13,-65638.2,-49.7858], Tmin=(1085.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-517.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsOsHH) + group(Cs-CsHHH) + polycyclic(s2_3_5_ane) + radical(CCsJOOC) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3801.71,'J/mol'), sigma=(6.35281,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.82 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.212558,0.0946781,-9.55049e-05,5.18623e-08,-1.21422e-11,-53122.3,28.8909], Tmin=(100,'K'), Tmax=(985.92,'K')), NASAPolynomial(coeffs=[11.2985,0.049701,-2.70755e-05,5.59098e-09,-4.09123e-13,-55308.3,-24.435], Tmin=(985.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-442.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'F[C](F)C1[C](F)C2OOC21F(1199)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {6,S} {7,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {8,S} {9,S}
11 C u1 p0 c0 {3,S} {4,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-420.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,190,488,555,1236,1407,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83322,0.0786775,-7.63173e-05,3.7489e-08,-7.7598e-12,-50469.4,24.4986], Tmin=(100,'K'), Tmax=(1112.64,'K')), NASAPolynomial(coeffs=[12.2317,0.0376989,-2.10716e-05,4.3868e-09,-3.21958e-13,-53005.8,-31.709], Tmin=(1112.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'FC1(F)CC2(F)C3OOC3(F)C12(1178)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {7,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {5,S} {7,S} {10,S} {14,S}
10 C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
11 C u0 p0 c0 {7,S} {12,S} {15,S} {16,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-713.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519147,0.0858849,-6.78272e-05,2.49904e-08,-3.78517e-12,-85717.9,17.8954], Tmin=(100,'K'), Tmax=(1470.06,'K')), NASAPolynomial(coeffs=[16.6,0.0421299,-2.31818e-05,4.74418e-09,-3.42127e-13,-90445.9,-65.8817], Tmin=(1470.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-713.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) - ring(Cyclobutane)"""),
)

species(
    label = 'CC(F)(F)C1=C(F)C2OOC12F(1200)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
10 C u0 p0 c0 {9,S} {14,S} {15,S} {16,S}
11 C u0 p0 c0 {7,S} {9,S} {12,D}
12 C u0 p0 c0 {4,S} {8,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-737.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.100685,0.0965929,-9.13941e-05,4.23999e-08,-8.06079e-12,-88604.9,25.2435], Tmin=(100,'K'), Tmax=(1226,'K')), NASAPolynomial(coeffs=[17.3735,0.0395814,-2.16418e-05,4.47084e-09,-3.26563e-13,-92889.6,-62.6201], Tmin=(1226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'CC(F)(F)C1C(F)=C2OOC21F(1016)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
10 C u0 p0 c0 {9,S} {14,S} {15,S} {16,S}
11 C u0 p0 c0 {6,S} {8,S} {12,D}
12 C u0 p0 c0 {4,S} {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-684.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.090024,0.101443,-0.000124361,8.86944e-08,-2.73254e-11,-82188.7,25.4501], Tmin=(100,'K'), Tmax=(769.011,'K')), NASAPolynomial(coeffs=[9.3143,0.0525272,-2.89514e-05,5.98398e-09,-4.37411e-13,-83635.1,-17.4505], Tmin=(769.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane)"""),
)

species(
    label = '[CH2]C(F)(F)C1[C](F)OOC=C1F(1201)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u1 p0 c0 {4,S} {5,S} {7,S}
11 C u0 p0 c0 {6,S} {9,D} {14,S}
12 C u1 p0 c0 {8,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-527.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,215,315,519,588,595,1205,1248,323,467,575,827,1418,395,473,707,1436,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744264,0.109684,-0.000138449,9.26785e-08,-2.49263e-11,-63241.1,32.2319], Tmin=(100,'K'), Tmax=(904.815,'K')), NASAPolynomial(coeffs=[15.8062,0.0365168,-1.71512e-05,3.30595e-09,-2.32475e-13,-66236.1,-45.9592], Tmin=(904.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(CsCsF1sO2s) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[CH2]C(F)(F)C=C(F)C1OO[C]1F(1202)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 C u1 p0 c0 {4,S} {6,S} {7,S}
12 C u1 p0 c0 {8,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-460.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,248,333,466,604,684,796,1061,1199,323,467,575,827,1418,3010,987.5,1337.5,450,1655,395,473,707,1436,3000,3100,440,815,1455,1000,180,180,180,432.983,729.391,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.171118,'amu*angstrom^2'), symmetry=1, barrier=(3.93434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171118,'amu*angstrom^2'), symmetry=1, barrier=(3.93434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171118,'amu*angstrom^2'), symmetry=1, barrier=(3.93434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.678068,0.112464,-0.000158364,1.21644e-07,-3.81229e-11,-55174.7,39.1158], Tmin=(100,'K'), Tmax=(775.026,'K')), NASAPolynomial(coeffs=[12.9599,0.0420753,-2.21275e-05,4.45244e-09,-3.193e-13,-57288.6,-23.2033], Tmin=(775.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-460.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFF) + group(CsCFHO) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1C(F)=CC1(F)O[O](1203)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {14,S}
12 C u1 p0 c0 {9,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-553.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,323,467,575,827,1418,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.86898,0.112879,-0.00013961,8.84605e-08,-2.24329e-11,-66452,34.0302], Tmin=(100,'K'), Tmax=(957.369,'K')), NASAPolynomial(coeffs=[17.7652,0.0350249,-1.76308e-05,3.52218e-09,-2.53271e-13,-70020,-55.0577], Tmin=(957.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-553.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C(F)(F)C1=C(F)C2OOC12F(1204)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {6,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {10,S} {12,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 C u1 p0 c0 {9,S} {14,S} {15,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-526.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,248,333,466,604,684,796,1061,1199,323,467,575,827,1418,3000,3100,440,815,1455,1000,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (183.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.131435,0.100809,-0.000109482,6.01958e-08,-1.36749e-11,-63182.9,26.6744], Tmin=(100,'K'), Tmax=(1038.49,'K')), NASAPolynomial(coeffs=[15.3529,0.0411655,-2.33313e-05,4.88978e-09,-3.60537e-13,-66398.9,-48.6133], Tmin=(1038.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-526.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1C(F)=C2OOC21F(1205)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
10 C u0 p0 c0 {6,S} {8,S} {11,D}
11 C u0 p0 c0 {4,S} {7,S} {10,D}
12 C u1 p0 c0 {9,S} {14,S} {15,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-472.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,323,467,575,827,1418,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (183.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.265364,0.1052,-0.000140013,1.04925e-07,-3.31301e-11,-56676.7,27.709], Tmin=(100,'K'), Tmax=(759.03,'K')), NASAPolynomial(coeffs=[10.5956,0.0479634,-2.69012e-05,5.57671e-09,-4.07565e-13,-58325.5,-21.6943], Tmin=(759.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    label = 'C=C(F)C1[C](F)C2OOC21F(1206)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {13,S}
9  C u1 p0 c0 {2,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {10,D} {14,S} {15,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-312.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,212,367,445,1450,323,467,575,827,1418,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.351146,0.0803986,-6.50329e-05,2.48278e-08,-3.81718e-12,-37465.1,25.8746], Tmin=(100,'K'), Tmax=(1500.75,'K')), NASAPolynomial(coeffs=[18.501,0.0320233,-1.66821e-05,3.34938e-09,-2.39242e-13,-42912.8,-69.0563], Tmin=(1500.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CdCsCdF) + group(Cds-CdsHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s)"""),
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
    label = '[CH2]C(F)(F)C1=C(F)C2OO[C]12(1207)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {4,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {6,S} {8,D}
11 C u1 p0 c0 {7,S} {13,S} {14,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-95.2311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,248,333,466,604,684,796,1061,1199,323,467,575,827,1418,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64465,0.0827614,-8.01672e-05,3.92322e-08,-8.05709e-12,-11340.4,27.856], Tmin=(100,'K'), Tmax=(1124.03,'K')), NASAPolynomial(coeffs=[12.9809,0.0388611,-2.1583e-05,4.48567e-09,-3.28954e-13,-14113.6,-33.1017], Tmin=(1124.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.2311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOO) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)=C1[C](F)C2OOC12F(1208)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {8,D} {11,S}
11 C u1 p0 c0 {10,S} {13,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-226.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,346,659,817,1284,271,519,563,612,1379,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.168092,0.0840177,-7.31112e-05,2.99001e-08,-4.89562e-12,-27113.3,24.0714], Tmin=(100,'K'), Tmax=(1428.12,'K')), NASAPolynomial(coeffs=[19.5418,0.029754,-1.61162e-05,3.29397e-09,-2.38066e-13,-32646.9,-76.2996], Tmin=(1428.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ane) + radical(CsCdCsF1s) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1[C](F)C2=C1OO2(1209)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u1 p0 c0 {3,S} {6,S} {10,S}
10 C u0 p0 c0 {5,S} {8,D} {9,S}
11 C u1 p0 c0 {7,S} {13,S} {14,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-0.59807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,346,659,817,1284,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.127749,0.079023,-6.72033e-05,2.65905e-08,-4.14644e-12,72.12,30.8578], Tmin=(100,'K'), Tmax=(1519.28,'K')), NASAPolynomial(coeffs=[21.2788,0.023336,-1.22229e-05,2.46482e-09,-1.76517e-13,-6354.75,-80.0298], Tmin=(1519.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.59807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + polycyclic(s2_4_4_ene_m) + radical(CsCdCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)[C]1C(F)C2OOC12F(1210)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-464.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0371676,0.0890763,-7.35567e-05,2.84304e-08,-4.41004e-12,-55744.2,30.787], Tmin=(100,'K'), Tmax=(1492.49,'K')), NASAPolynomial(coeffs=[20.5334,0.0339453,-1.81483e-05,3.68047e-09,-2.64291e-13,-61884.4,-76.6915], Tmin=(1492.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-464.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1C(F)[C]2OOC21F(1211)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {7,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {6,S} {8,S} {9,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-468.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,316,385,515,654,689,1295,259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0572858,0.100068,-0.0001031,5.5547e-08,-1.26363e-11,-56266.8,26.1534], Tmin=(100,'K'), Tmax=(1024.34,'K')), NASAPolynomial(coeffs=[13.4204,0.0474389,-2.6032e-05,5.38912e-09,-3.94819e-13,-59028,-39.1927], Tmin=(1024.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'CC(F)(F)[C]1[C](F)C2OOC12F(1212)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
10 C u0 p0 c0 {9,S} {14,S} {15,S} {16,S}
11 C u1 p0 c0 {7,S} {9,S} {12,S}
12 C u1 p0 c0 {4,S} {8,S} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-486.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.669876,0.0838251,-6.78702e-05,2.66657e-08,-4.41415e-12,-58375,29.3115], Tmin=(100,'K'), Tmax=(1337.21,'K')), NASAPolynomial(coeffs=[13.6309,0.0450547,-2.43799e-05,4.98355e-09,-3.60542e-13,-61841.3,-36.9842], Tmin=(1337.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-486.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'CC(F)(F)C1[C](F)[C]2OOC21F(1213)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
10 C u0 p0 c0 {9,S} {14,S} {15,S} {16,S}
11 C u1 p0 c0 {6,S} {8,S} {12,S}
12 C u1 p0 c0 {4,S} {7,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-490.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.247025,0.105134,-0.000132399,9.8127e-08,-3.14358e-11,-58858.5,27.9117], Tmin=(100,'K'), Tmax=(741.528,'K')), NASAPolynomial(coeffs=[9.11932,0.0546091,-3.01921e-05,6.23753e-09,-4.55524e-13,-60247.6,-14.4744], Tmin=(741.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-490.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = '[CH2]C(F)(F)C1[C]2OOC2C1(F)F(1214)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-457.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,222,329,445,522,589,1214,1475,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.119658,0.101247,-0.000110524,6.47878e-08,-1.60466e-11,-54834.8,28.4378], Tmin=(100,'K'), Tmax=(949.479,'K')), NASAPolynomial(coeffs=[12.605,0.0476397,-2.58334e-05,5.32228e-09,-3.89003e-13,-57251.1,-32.2912], Tmin=(949.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-457.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'FCC(F)(F)C1[C](F)C2OO[C]21(1215)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {5,S} {11,S} {12,S} {14,S}
10 C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
11 C u1 p0 c0 {6,S} {7,S} {9,S}
12 C u1 p0 c0 {4,S} {7,S} {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-415.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,222,329,445,522,589,1214,1475,528,1116,1182,1331,1402,1494,3075,3110,212,367,445,1450,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0192414,0.0994153,-0.00011621,7.88554e-08,-2.31879e-11,-49857,30.1064], Tmin=(100,'K'), Tmax=(802.153,'K')), NASAPolynomial(coeffs=[9.40596,0.052416,-2.8323e-05,5.81337e-09,-4.2361e-13,-51369.1,-13.287], Tmin=(802.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-415.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = '[CH2][C](F)C1C(F)(F)C2OOC21F(1216)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {10,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
11 C u1 p0 c0 {4,S} {7,S} {12,S}
12 C u1 p0 c0 {11,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-470.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.134084,0.102214,-0.000111078,6.42497e-08,-1.57101e-11,-56443,28.1113], Tmin=(100,'K'), Tmax=(959.6,'K')), NASAPolynomial(coeffs=[12.8451,0.0481109,-2.65064e-05,5.49482e-09,-4.02904e-13,-58933.9,-33.9702], Tmin=(959.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-470.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = 'FC[C](F)C1[C](F)C2OOC21F(1082)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {2,S} {12,S} {15,S} {16,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 C u1 p0 c0 {3,S} {7,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-428.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,551,1088,1226,1380,1420,1481,3057,3119,165,259,333,401,397,493,1395,1505,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.11947,0.0970573,-0.000102499,5.93615e-08,-1.48558e-11,-51456.4,29.361], Tmin=(100,'K'), Tmax=(928.357,'K')), NASAPolynomial(coeffs=[10.8398,0.0508659,-2.78636e-05,5.7642e-09,-4.22231e-13,-53446.8,-21.5615], Tmin=(928.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-428.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(CsCsCsF1s)"""),
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
    E0 = (-151.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (154.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (7.98515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (19.1177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (7.98515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (116.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (274.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-143.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-88.5498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-126.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-83.1727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-21.7437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-151.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-37.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (10.2256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (60.9746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (125.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-20.8025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (228.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (157.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (53.8874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (213.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (22.1047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-31.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-57.9456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-108.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (79.8026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (115.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (73.0868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (20.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['CH2CF2(57)', 'FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C1[C](F)C(F)C12OO2(1196)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['[CH2]C(F)(F)[CH]C1(F)C2OOC21F(1197)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['[CH2]C(F)(F)C1C2(F)[CH]OOC12F(1198)'],
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
    reactants = ['CH2(T)(18)', 'F[C](F)C1[C](F)C2OOC21F(1199)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['FC1(F)CC2(F)C3OOC3(F)C12(1178)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['CC(F)(F)C1=C(F)C2OOC12F(1200)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['CC(F)(F)C1C(F)=C2OOC21F(1016)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(F)(F)C1[C](F)OOC=C1F(1201)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c6_gamma;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(F)(F)C=C(F)C1OO[C]1F(1202)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(F)(F)C1C(F)=CC1(F)O[O](1203)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.10864e+11,'s^-1'), n=0.310877, Ea=(88.8126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 88.0 to 88.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](F)F(163)', 'FC1=CC2(F)OOC12(939)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.05802e-06,'m^3/(mol*s)'), n=3.12542, Ea=(6.81068,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3001154165706074, var=1.9196358956114195, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-2R!H-R_Ext-2R!H-R_Ext-3C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-2R!H-R_Ext-2R!H-R_Ext-3C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[CH2]C(F)(F)C1=C(F)C2OOC12F(1204)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(205.9,'m^3/(mol*s)'), n=1.596, Ea=(11.7533,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[CH2]C(F)(F)C1C(F)=C2OOC21F(1205)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(45.45,'m^3/(mol*s)'), n=1.71, Ea=(8.43039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_N-Sp-9R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_N-Sp-9R!H=5R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'C=C(F)C1[C](F)C2OOC21F(1206)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(52.1927,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2CF2(57)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(3.568,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][C](F)F(163)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', '[CH2]C(F)(F)C1=C(F)C2OO[C]12(1207)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(220.965,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[CH2]C(F)=C1[C](F)C2OOC12F(1208)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(248.409,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[CH2]C(F)(F)C1[C](F)C2=C1OO2(1209)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(181.713,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['[CH2]C(F)(F)[C]1C(F)C2OOC12F(1210)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.93363e+09,'s^-1'), n=1.033, Ea=(174.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_Cs2] for rate rule [R2H_S_cy4;C_rad_out_noH;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(F)(F)C1C(F)[C]2OOC21F(1211)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.26888e+09,'s^-1'), n=0.972714, Ea=(123.905,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_NDMustO;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['CC(F)(F)[C]1[C](F)C2OOC12F(1212)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.4313e+06,'s^-1'), n=1.72764, Ea=(94.0044,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['CC(F)(F)C1[C](F)[C]2OOC21F(1213)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(208392,'s^-1'), n=1.49392, Ea=(43.3506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_CCC;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(F)(F)C1[C]2OOC2C1(F)F(1214)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(223.716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction28',
    reactants = ['FCC(F)(F)C1[C](F)C2OO[C]21(1215)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(218.416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    products = ['[CH2][C](F)C1C(F)(F)C2OOC21F(1216)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(225.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['FC[C](F)C1[C](F)C2OOC21F(1082)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(136.48,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #214',
    isomers = [
        '[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)',
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
    label = 'PDepNetwork #214',
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

