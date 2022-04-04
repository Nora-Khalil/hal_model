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
    label = '[CH2]C(F)(F)C12[CH]C(F)(OO1)C2F(1217)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {7,S} {9,S} {14,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-502.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.796984,0.101806,-9.78327e-05,4.48821e-08,-8.14781e-12,-60312.5,29.9271], Tmin=(100,'K'), Tmax=(1318.75,'K')), NASAPolynomial(coeffs=[23.1704,0.0291092,-1.5145e-05,3.08135e-09,-2.23557e-13,-66634,-92.3334], Tmin=(1318.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-502.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsOs) + group(CsCCFO) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s3_4_5_ane) + radical(CCJCOOH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)[CH]C2(OO2)C1F(1218)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {7,S} {8,S} {14,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-443.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,250,417,511,1155,1315,1456,3119,215,315,519,588,595,1205,1248,2950,1000,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80847,0.117564,-0.000130116,6.80638e-08,-1.37053e-11,-53178.7,35.7245], Tmin=(100,'K'), Tmax=(1222.21,'K')), NASAPolynomial(coeffs=[28.5084,0.018344,-8.34391e-06,1.64138e-09,-1.18664e-13,-60589.4,-116.62], Tmin=(1222.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s1_3_4_ane) + radical(CCJCOOH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)(F)[C](F)C1C2OOC21F(1219)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
10 C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
11 C u1 p0 c0 {4,S} {7,S} {10,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-453.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.424162,0.105091,-0.00011315,6.17368e-08,-1.37913e-11,-54345.5,30.7794], Tmin=(100,'K'), Tmax=(1062.41,'K')), NASAPolynomial(coeffs=[16.7205,0.0405395,-2.20083e-05,4.54328e-09,-3.32554e-13,-57988.3,-52.9713], Tmin=(1062.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-453.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_3_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,215,315,519,588,595,1205,1248,212,367,445,1450,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3827.89,'J/mol'), sigma=(6.67154,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=597.91 K, Pc=29.25 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.141016,0.100934,-0.000104417,5.5609e-08,-1.23631e-11,-55799.9,28.3539], Tmin=(100,'K'), Tmax=(1053.25,'K')), NASAPolynomial(coeffs=[14.6145,0.0448961,-2.46116e-05,5.09555e-09,-3.73346e-13,-58908.2,-43.5986], Tmin=(1053.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-465.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)C2OO[C](F)C21(1220)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {5,S} {7,S} {8,S} {14,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {4,S} {6,S} {8,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-521.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22954,0.101468,-0.000107202,5.8376e-08,-1.31556e-11,-62616.1,29.3057], Tmin=(100,'K'), Tmax=(1047.24,'K')), NASAPolynomial(coeffs=[15.1815,0.042606,-2.28927e-05,4.70649e-09,-3.43724e-13,-65844,-45.7551], Tmin=(1047.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-521.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(Cs-CsHHH) + polycyclic(s2_3_5_ane) + radical(CsCsF1sO2s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2750,3150,900,1100,316,385,515,654,689,1295,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3875.66,'J/mol'), sigma=(6.72331,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.37 K, Pc=28.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0223281,0.100306,-0.000109539,6.5089e-08,-1.65562e-11,-53853.1,29.1069], Tmin=(100,'K'), Tmax=(919.545,'K')), NASAPolynomial(coeffs=[11.5124,0.0501287,-2.76855e-05,5.74367e-09,-4.21305e-13,-55974.3,-25.5738], Tmin=(919.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-448.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFFH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'F[C](F)C1(F)[CH]C2(F)OOC21(1221)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {13,S}
11 C u1 p0 c0 {3,S} {4,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-404.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,289,311,382,485,703,1397,316,385,515,654,689,1295,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481404,0.0879152,-0.000104139,6.70662e-08,-1.83371e-11,-48538.8,25.0537], Tmin=(100,'K'), Tmax=(864.084,'K')), NASAPolynomial(coeffs=[10.5022,0.0415287,-2.36169e-05,4.94306e-09,-3.63983e-13,-50270.6,-21.8271], Tmin=(864.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(Csj(Cs-F1sCsCs)(F1s)(F1s))"""),
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
    label = '[CH2]C(F)(F)C1(F)[CH]OOC(F)=C1(1222)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {5,S} {7,S} {14,S}
11 C u0 p0 c0 {4,S} {6,S} {9,D}
12 C u1 p0 c0 {8,S} {15,S} {16,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-534.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2750,3150,900,1100,326,540,652,719,1357,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.9892,0.118544,-0.00017756,1.473e-07,-4.91451e-11,-64166.8,32.23], Tmin=(100,'K'), Tmax=(790.599,'K')), NASAPolynomial(coeffs=[12.5273,0.0431077,-2.10593e-05,4.05285e-09,-2.81386e-13,-66083.7,-28.4101], Tmin=(790.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-534.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFO) + ring(34dihydro12dioxin) + radical(CCsJOOC) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C(F)=CC1(F)[CH]OO1(1223)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {10,S} {12,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {4,S} {8,S} {9,D}
11 C u1 p0 c0 {6,S} {7,S} {14,S}
12 C u1 p0 c0 {8,S} {15,S} {16,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-469.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,248,333,466,604,684,796,1061,1199,3010,987.5,1337.5,450,1655,323,467,575,827,1418,2950,1000,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1462,0.123818,-0.000195438,1.68751e-07,-5.85041e-11,-56290.8,36.8741], Tmin=(100,'K'), Tmax=(770.375,'K')), NASAPolynomial(coeffs=[12.4255,0.0443681,-2.32499e-05,4.60743e-09,-3.25042e-13,-58115.3,-23.3308], Tmin=(770.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCsJOO) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)C=C(F)C1O[O](1224)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 C u1 p0 c0 {9,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-519.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,3150,900,1100,215,315,519,588,595,1205,1248,323,467,575,827,1418,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.891622,0.114107,-0.000146566,9.69235e-08,-2.56317e-11,-62366.7,36.0272], Tmin=(100,'K'), Tmax=(920.078,'K')), NASAPolynomial(coeffs=[17.2953,0.0350396,-1.76613e-05,3.5213e-09,-2.52574e-13,-65713.3,-50.1992], Tmin=(920.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-519.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C(F)(F)C1=CC2(F)OOC12(1225)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {11,S}
9  C u0 p0 c0 {6,S} {8,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {13,S}
11 C u1 p0 c0 {8,S} {14,S} {15,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-318.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,1380,1390,370,380,2900,435,248,333,466,604,684,796,1061,1199,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288431,0.0878587,-8.0407e-05,3.56791e-08,-6.49236e-12,-38118.8,25.043], Tmin=(100,'K'), Tmax=(1271.81,'K')), NASAPolynomial(coeffs=[16.4688,0.036969,-2.03864e-05,4.21693e-09,-3.07805e-13,-42234.5,-56.9086], Tmin=(1271.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
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
    label = '[CH2]C(F)(F)C1(F)C=C2OOC21(1226)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {4,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
9  C u0 p0 c0 {5,S} {7,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {13,S}
11 C u1 p0 c0 {8,S} {14,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-267.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320516,0.0874075,-8.59212e-05,4.25605e-08,-8.6991e-12,-32095.2,27.1646], Tmin=(100,'K'), Tmax=(1145.29,'K')), NASAPolynomial(coeffs=[14.7756,0.0369212,-1.97976e-05,4.0697e-09,-2.96994e-13,-35406.3,-44.5337], Tmin=(1145.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'C=C(F)C1(F)[CH]C2(F)OOC21(1227)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {6,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
9  C u1 p0 c0 {6,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {10,D} {14,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-320.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,316,385,515,654,689,1295,323,467,575,827,1418,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533564,0.0844006,-7.57016e-05,3.38419e-08,-6.33232e-12,-38481.6,23.8905], Tmin=(100,'K'), Tmax=(1221.2,'K')), NASAPolynomial(coeffs=[13.8866,0.0406634,-2.19795e-05,4.51462e-09,-3.28566e-13,-41743,-43.1988], Tmin=(1221.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CdCsCdF) + group(Cds-CdsHH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH)"""),
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
    label = '[CH2]C(F)=C1[CH]C2(F)OOC12(1228)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,D}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {7,D} {10,S}
10 C u1 p0 c0 {9,S} {13,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-46.6775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,316,385,515,654,689,1295,271,519,563,612,1379,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.136254,0.0745271,-3.94226e-05,-6.69464e-09,7.82749e-12,-5451.13,24.0891], Tmin=(100,'K'), Tmax=(1097.26,'K')), NASAPolynomial(coeffs=[22.6538,0.0227977,-1.15634e-05,2.41743e-09,-1.81205e-13,-12339.7,-96.5747], Tmin=(1097.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.6775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ane) + radical(C=CCJCO) + radical(Csj(Cd-F1sCd)(H)(H))"""),
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
    label = '[CH2]C(F)(F)C1=CC2(F)OO[C]12(1229)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {5,S} {6,S} {8,S}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 C u1 p0 c0 {7,S} {13,S} {14,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-131.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,248,333,466,604,684,796,1061,1199,2950,1000,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.393912,0.0907786,-0.000104356,6.60885e-08,-1.80101e-11,-15681.1,25.5069], Tmin=(100,'K'), Tmax=(860.79,'K')), NASAPolynomial(coeffs=[10.0612,0.0458562,-2.60756e-05,5.46219e-09,-4.02515e-13,-17345.5,-19.6829], Tmin=(860.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOO) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)[CH]C2=C1OO2(1230)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {12,S}
10 C u0 p0 c0 {5,S} {8,D} {9,S}
11 C u1 p0 c0 {7,S} {13,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (50.8074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461762,0.0784532,-7.0052e-05,3.00221e-08,-5.18147e-12,6237.12,31.0132], Tmin=(100,'K'), Tmax=(1357.87,'K')), NASAPolynomial(coeffs=[17.3162,0.0288037,-1.52058e-05,3.09457e-09,-2.23806e-13,1659.89,-55.4557], Tmin=(1357.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.8074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + polycyclic(s2_4_4_ene_m) + radical(CCJCO) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)CC2(F)OO[C]21(1231)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {11,S}
9  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-465.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.241175,0.105024,-0.000119012,7.30024e-08,-1.89716e-11,-55803.1,27.0574], Tmin=(100,'K'), Tmax=(906.995,'K')), NASAPolynomial(coeffs=[12.3156,0.0496465,-2.74279e-05,5.68584e-09,-4.16787e-13,-58080.9,-32.2962], Tmin=(906.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-465.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'CC(F)(F)C1(F)[CH]C2(F)OO[C]21(1232)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
10 C u0 p0 c0 {8,S} {13,S} {14,S} {15,S}
11 C u1 p0 c0 {6,S} {7,S} {9,S}
12 C u1 p0 c0 {7,S} {9,S} {16,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-476.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27805,0.106137,-0.000129478,8.97797e-08,-2.66747e-11,-57203.7,27.8617], Tmin=(100,'K'), Tmax=(797.825,'K')), NASAPolynomial(coeffs=[10.3522,0.0528437,-2.92863e-05,6.06386e-09,-4.43703e-13,-58900.1,-21.0228], Tmin=(797.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-476.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(F)(F)[C]1C(F)C2(F)OOC12(1233)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {6,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-439.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2950,1000,259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.129177,0.0887224,-7.29388e-05,2.81e-08,-4.32496e-12,-52726.2,32.4297], Tmin=(100,'K'), Tmax=(1512.84,'K')), NASAPolynomial(coeffs=[21.2522,0.0321893,-1.68854e-05,3.39874e-09,-2.43019e-13,-59195.5,-79.5745], Tmin=(1512.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-439.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCC(F)(F)[C]1[CH]C2(F)OOC12(1234)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
10 C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
11 C u1 p0 c0 {7,S} {9,S} {12,S}
12 C u1 p0 c0 {8,S} {11,S} {16,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-432.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,316,385,515,654,689,1295,215,315,519,588,595,1205,1248,528,1116,1182,1331,1402,1494,3075,3110,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0587558,0.0852384,-6.66155e-05,2.40724e-08,-3.47066e-12,-51935,33.3177], Tmin=(100,'K'), Tmax=(1599.92,'K')), NASAPolynomial(coeffs=[21.2102,0.0323573,-1.70369e-05,3.41358e-09,-2.42562e-13,-58703.1,-78.6658], Tmin=(1599.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-432.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)C(F)[C]2OOC21(1235)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
11 C u1 p0 c0 {6,S} {8,S} {9,S}
12 C u1 p0 c0 {10,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-408.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([321,327,383,445,494,1500,2950,1000,259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.246992,0.104125,-0.00012092,7.69607e-08,-2.0662e-11,-48993.3,29.1914], Tmin=(100,'K'), Tmax=(883.131,'K')), NASAPolynomial(coeffs=[12.0835,0.048276,-2.60593e-05,5.3505e-09,-3.90112e-13,-51171.2,-28.7635], Tmin=(883.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-408.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCC(F)(F)C1(F)[CH][C]2OOC21(1236)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
10 C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
11 C u1 p0 c0 {6,S} {8,S} {12,S}
12 C u1 p0 c0 {7,S} {11,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-401.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2750,3150,900,1100,222,329,445,522,589,1214,1475,528,1116,1182,1331,1402,1494,3075,3110,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0857052,0.100924,-0.000115516,7.41529e-08,-2.03809e-11,-48200.8,30.1779], Tmin=(100,'K'), Tmax=(857.977,'K')), NASAPolynomial(coeffs=[10.7368,0.0504688,-2.73057e-05,5.61221e-09,-4.09434e-13,-50057.9,-20.3765], Tmin=(857.977,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2][C](F)C1(F)C(F)C2(F)OOC21(1237)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {7,S} {9,S} {14,S}
11 C u1 p0 c0 {4,S} {7,S} {12,S}
12 C u1 p0 c0 {11,S} {15,S} {16,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-409.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2950,1000,333,384,448,608,1254,1480,250,417,511,1155,1315,1456,3119,212,367,445,1450,3000,3100,440,815,1455,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.253617,0.104733,-0.000121543,7.70895e-08,-2.06578e-11,-49099.2,29.8627], Tmin=(100,'K'), Tmax=(883.533,'K')), NASAPolynomial(coeffs=[12.0933,0.0488385,-2.66552e-05,5.49711e-09,-4.01757e-13,-51281.1,-28.1762], Tmin=(883.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-409.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC[C](F)C1(F)[CH]C2(F)OOC21(1120)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {2,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {12,S} {14,S} {15,S}
11 C u1 p0 c0 {7,S} {9,S} {16,S}
12 C u1 p0 c0 {4,S} {7,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-427.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2750,3150,900,1100,316,385,515,654,689,1295,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (184.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0432941,0.100801,-0.00011005,6.52317e-08,-1.65361e-11,-51307,29.4965], Tmin=(100,'K'), Tmax=(922.799,'K')), NASAPolynomial(coeffs=[11.6414,0.0501519,-2.77186e-05,5.75181e-09,-4.21945e-13,-53463.5,-25.9362], Tmin=(922.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-427.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(CsCsCsF1s)"""),
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
    E0 = (-170.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (128.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (136.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-12.7436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-12.7436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-12.7436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-12.7436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (258.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-161.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-122.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-61.9597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-170.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (93.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-63.2573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (130.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (91.4203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-50.7527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (196.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (246.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (112.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (218.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-21.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-75.9797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (25.5871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (68.5329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (40.5221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (54.1807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (79.1877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-4.08538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['CH2CF2(57)', 'FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['[CH2]C(F)(F)C12[CH]C(F)(OO1)C2F(1217)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(OO2)C1F(1218)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['[CH2]C(F)(F)[C](F)C1C2OOC21F(1219)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['[CH2]C(F)(F)C1[C](F)C2OOC21F(1012)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['[CH2]C(F)(F)C1(F)C2OO[C](F)C21(1220)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['F[C](F)CC1(F)[CH]C2(F)OOC12(1013)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(18)', 'F[C](F)C1(F)[CH]C2(F)OOC21(1221)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['FC12OOC1C1(F)C2CC1(F)F(1135)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]OOC(F)=C1(1222)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra_pri;radadd_intra_csHNd] for rate rule [Rn0c6_gamma;doublebond_intra_pri;radadd_intra_csHO]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(F)(F)C(F)=CC1(F)[CH]OO1(1223)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.23261e+08,'s^-1'), n=1.06333, Ea=(126.217,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_csHNd] for rate rule [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(F)(F)C1(F)C=C(F)C1O[O](1224)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.83633e+10,'s^-1'), n=0.345735, Ea=(68.5957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 68.2 to 68.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[CH2]C(F)(F)C1=CC2(F)OOC12(1225)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(57.1836,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](F)F(163)', 'FC1=CC2(F)OOC12(939)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.000202,'m^3/(mol*s)'), n=2.41, Ea=(13.0857,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-2R!H-R_Ext-1R!H-R_Ext-3C-R_Ext-3C-R_Sp-6R!H-2R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-2R!H-R_Ext-1R!H-R_Ext-3C-R_Ext-3C-R_Sp-6R!H-2R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[CH2]C(F)(F)C1(F)C=C2OOC21(1226)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(44.0787,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'C=C(F)C1(F)[CH]C2(F)OOC21(1227)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(58.1543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2CF2(57)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(5.47907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C](F)F(163)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F2(78)', '[CH2]C(F)=C1[CH]C2(F)OOC12(1228)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(20.3085,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', '[CH2]C(F)(F)C1=CC2(F)OO[C]12(1229)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(243.699,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[CH2]C(F)(F)C1(F)[CH]C2=C1OO2(1230)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(167.19,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['[CH2]C(F)(F)C1(F)CC2(F)OO[C]21(1231)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.05815e+09,'s^-1'), n=0.95, Ea=(148.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    products = ['CC(F)(F)C1(F)[CH]C2(F)OO[C]21(1232)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.74067e+06,'s^-1'), n=1.70063, Ea=(94.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cs_H_out_OOH/Cs] + [R4H_SSS;C_rad_out_2H;Cs_H_out] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(F)(F)[C]1C(F)C2(F)OOC12(1233)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(183.914,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FCC(F)(F)[C]1[CH]C2(F)OOC12(1234)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(220.219,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(F)(F)C1(F)C(F)[C]2OOC21(1235)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(167.765,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['FCC(F)(F)C1(F)[CH][C]2OOC21(1236)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(174.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C](F)C1(F)C(F)C2(F)OOC21(1237)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(207.309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['FC[C](F)C1(F)[CH]C2(F)OOC21(1120)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(142.326,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #213',
    isomers = [
        '[CH2]C(F)(F)C1(F)[CH]C2(F)OOC21(1011)',
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
    label = 'PDepNetwork #213',
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

