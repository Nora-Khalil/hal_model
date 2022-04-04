species(
    label = '[O]C(F)C1(F)[CH]C2(F)OOC21(1023)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
11 C u1 p0 c0 {7,S} {9,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-366.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2750,3150,900,1100,316,385,515,654,689,1295,391,562,707,872,1109,1210,1289,3137,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.553865,0.0821399,-7.54375e-05,3.4009e-08,-6.32621e-12,-44012.8,25.6287], Tmin=(100,'K'), Tmax=(1242.14,'K')), NASAPolynomial(coeffs=[14.9127,0.0359008,-1.95995e-05,4.0403e-09,-2.94561e-13,-47579.9,-46.7581], Tmin=(1242.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCOOH)"""),
)

species(
    label = 'CHFO(47)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.294,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C12[CH]C(F)(OO1)C2F(1431)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
11 C u1 p0 c0 {7,S} {9,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-405.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.341046,0.0834029,-7.32733e-05,3.00159e-08,-4.7867e-12,-48655.9,28.5519], Tmin=(100,'K'), Tmax=(1511.79,'K')), NASAPolynomial(coeffs=[23.9622,0.019099,-9.47009e-06,1.87968e-09,-1.33846e-13,-56004.1,-98.7413], Tmin=(1511.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-405.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(CsCCFO) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFHO) + polycyclic(s3_4_5_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCOOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)C1(F)[CH]C2(OO2)C1F(1432)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {5,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
11 C u1 p0 c0 {7,S} {8,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-359.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,250,417,511,1155,1315,1456,3119,391,562,707,872,1109,1210,1289,3137,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.930858,0.095621,-9.5076e-05,4.13552e-08,-5.95561e-12,-43052,32.5321], Tmin=(100,'K'), Tmax=(1077.47,'K')), NASAPolynomial(coeffs=[25.5458,0.0145661,-6.23165e-06,1.23157e-09,-9.08845e-14,-49758.1,-101.821], Tmin=(1077.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-359.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s1_3_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCOOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)[C](F)C1C2OOC21F(1433)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {13,S}
9  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
10 C u0 p0 c0 {2,S} {6,S} {11,S} {14,S}
11 C u1 p0 c0 {3,S} {7,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-361.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,391,562,707,872,1109,1210,1289,3137,212,367,445,1450,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400149,0.0828659,-7.75946e-05,3.55913e-08,-6.65214e-12,-43333.9,27.9453], Tmin=(100,'K'), Tmax=(1252.12,'K')), NASAPolynomial(coeffs=[16.114,0.0326682,-1.74609e-05,3.57512e-09,-2.59911e-13,-47269.1,-51.3983], Tmin=(1252.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-361.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(F)C1[C](F)C2OOC21F(1024)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {6,S} {7,S} {14,S}
11 C u1 p0 c0 {3,S} {7,S} {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-368.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,391,562,707,872,1109,1210,1289,3137,212,367,445,1450,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4213.06,'J/mol'), sigma=(6.83843,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.07 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700842,0.0786221,-6.8804e-05,2.94181e-08,-5.20195e-12,-44161.7,25.5463], Tmin=(100,'K'), Tmax=(1294.88,'K')), NASAPolynomial(coeffs=[14.5168,0.0359432,-1.93642e-05,3.96409e-09,-2.8757e-13,-47739.7,-44.6779], Tmin=(1294.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(F)C1(F)C2OO[C](F)C21(1434)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {2,S} {6,S} {7,S} {14,S}
11 C u1 p0 c0 {3,S} {5,S} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-437.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660904,0.0795885,-7.31788e-05,3.37223e-08,-6.4568e-12,-52490.7,26.0503], Tmin=(100,'K'), Tmax=(1206.7,'K')), NASAPolynomial(coeffs=[13.7446,0.0362182,-1.92666e-05,3.93734e-09,-2.86019e-13,-55648.3,-39.5292], Tmin=(1206.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-437.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_5_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'F[CH]C1(F)[CH]C2(F)OOC12(1107)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u1 p0 c0 {3,S} {7,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-204.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,289,311,382,485,703,1397,316,385,515,654,689,1295,334,575,1197,1424,3202,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731605,0.0822555,-9.12762e-05,5.49849e-08,-1.42194e-11,-24478.9,21.9819], Tmin=(100,'K'), Tmax=(903.667,'K')), NASAPolynomial(coeffs=[9.97313,0.041348,-2.33725e-05,4.88909e-09,-3.60115e-13,-26149.2,-21.6668], Tmin=(903.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(Csj(Cs-F1sCsCs)(F1s)(H))"""),
)

species(
    label = 'FC1OC2C3(F)OOC3C12F(1028)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-578.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586035,0.074362,-5.47254e-05,1.7615e-08,-2.21491e-12,-69511.9,17.8853], Tmin=(100,'K'), Tmax=(1828.96,'K')), NASAPolynomial(coeffs=[23.0611,0.0252076,-1.44114e-05,2.92017e-09,-2.06246e-13,-77733,-104.113], Tmin=(1828.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-578.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCFHO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) - ring(Cyclobutane)"""),
)

species(
    label = 'O=C(F)C1(F)CC2(F)OOC21(1435)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
11 C u0 p0 c0 {3,S} {6,D} {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-750.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622673,0.0754073,-5.9427e-05,2.14685e-08,-3.10139e-12,-90171,22.9431], Tmin=(100,'K'), Tmax=(1586.81,'K')), NASAPolynomial(coeffs=[18.8047,0.0295741,-1.61009e-05,3.26583e-09,-2.33564e-13,-95941.3,-73.1699], Tmin=(1586.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-750.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(COCsFO) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = '[O]C(F)C1(F)[CH]OOC(F)=C1(1436)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {10,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {4,S} {7,S} {14,S}
11 C u0 p0 c0 {3,S} {5,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-450.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,2750,3150,900,1100,326,540,652,719,1357,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0130896,0.0954135,-0.00013817,1.145e-07,-3.8527e-11,-54044.5,28.6851], Tmin=(100,'K'), Tmax=(790.218,'K')), NASAPolynomial(coeffs=[9.94983,0.0387009,-1.8594e-05,3.56116e-09,-2.46904e-13,-55422.9,-15.7936], Tmin=(790.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFO) + ring(34dihydro12dioxin) + radical(O2sj(Cs-CsF1sH)) + radical(CCsJOOC)"""),
)

species(
    label = '[O]C(F)C(F)=CC1(F)[CH]OO1(1437)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {10,S} {12,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {3,S} {8,S} {9,D}
11 C u1 p0 c0 {5,S} {7,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-395.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,361,769,965,1078,1132,1246,3247,3010,987.5,1337.5,450,1655,323,467,575,827,1418,2950,1000,180,180,180,180,992.516,1600,1812.73,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.136632,'amu*angstrom^2'), symmetry=1, barrier=(3.14145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136632,'amu*angstrom^2'), symmetry=1, barrier=(3.14145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.110954,0.0927044,-0.000118315,8.19745e-08,-2.32599e-11,-47406.5,32.4332], Tmin=(100,'K'), Tmax=(851.415,'K')), NASAPolynomial(coeffs=[12.1958,0.0359283,-1.82875e-05,3.65115e-09,-2.61667e-13,-49464.3,-23.9254], Tmin=(851.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-395.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs-Cs(F)-O2s-O2s) + radical(O2sj(Cs-F1sCdH)) + radical(CCsJOO)"""),
)

species(
    label = '[O]OC1C(F)=CC1(F)C([O])F(1438)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {9,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {13,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u0 p0 c0 {3,S} {8,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-435.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,3150,900,1100,391,562,707,872,1109,1210,1289,3137,323,467,575,827,1418,296.867,296.867,296.867,296.867,296.867,296.867,2279.82,4000],'cm^-1')),
        HinderedRotor(inertia=(0.494788,'amu*angstrom^2'), symmetry=1, barrier=(30.9435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494788,'amu*angstrom^2'), symmetry=1, barrier=(30.9435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.041788,0.0915678,-0.000109664,6.79711e-08,-1.69424e-11,-52242.6,32.6297], Tmin=(100,'K'), Tmax=(971.309,'K')), NASAPolynomial(coeffs=[14.9341,0.030241,-1.49596e-05,2.97193e-09,-2.13189e-13,-55135.7,-38.7839], Tmin=(971.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sH)) + radical(ROOJ)"""),
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
    label = '[O]C(F)C1=CC2(F)OOC12(1439)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {5,S} {9,S} {12,S}
9  C u0 p0 c0 {6,S} {8,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-237.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,1380,1390,370,380,2900,435,361,769,965,1078,1132,1246,3247,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231302,0.073188,-6.00087e-05,2.23542e-08,-3.2414e-12,-28418.5,25.0111], Tmin=(100,'K'), Tmax=(1642.98,'K')), NASAPolynomial(coeffs=[22.8248,0.018182,-9.78958e-06,1.97698e-09,-1.40743e-13,-35842.6,-95.2071], Tmin=(1642.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CsCFHO) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(O2sj(Cs-F1sCdH))"""),
)

species(
    label = '[O][CH]F(388)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53591,0.00799828,-5.22385e-07,-4.56598e-09,2.08746e-12,-2348.33,9.31393], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97189,0.00388217,-1.62977e-06,3.36434e-10,-2.54185e-14,-3160.19,-3.94933], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O]C(F)C1(F)C=C2OOC21(1440)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {5,S} {6,S} {12,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-183.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,391,562,707,872,1109,1210,1289,3137,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1098,0.0666747,-5.56189e-05,2.22456e-08,-3.63164e-12,-21965.5,24.2738], Tmin=(100,'K'), Tmax=(1404.01,'K')), NASAPolynomial(coeffs=[14.3519,0.0289484,-1.53133e-05,3.10732e-09,-2.2385e-13,-25683.9,-44.105], Tmin=(1404.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'O=CC1(F)[CH]C2(F)OOC12(1441)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u0 p0 c0 {5,D} {7,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-296.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,1380,1390,370,380,2900,435,316,385,515,654,689,1295,2782.5,750,1395,475,1775,1000,499.806,499.81,499.81,499.811,499.812,499.817,499.818,499.818,499.822,499.831],'cm^-1')),
        HinderedRotor(inertia=(0.000674785,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06769,0.0666601,-5.0207e-05,1.70369e-08,-2.30286e-12,-35535.4,22.7232], Tmin=(100,'K'), Tmax=(1679.24,'K')), NASAPolynomial(coeffs=[17.6824,0.0270838,-1.48555e-05,3.00236e-09,-2.13469e-13,-41115.5,-66.0456], Tmin=(1679.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-296.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH)"""),
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
    label = 'O=C(F)C1(F)[CH]C2(F)OOC21(1442)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {6,D} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-550.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,1380,1390,370,380,2900,435,316,385,515,654,689,1295,486,617,768,1157,1926,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.962118,0.0737974,-6.30795e-05,2.55607e-08,-4.26323e-12,-66083.4,24.0278], Tmin=(100,'K'), Tmax=(1360.36,'K')), NASAPolynomial(coeffs=[14.5867,0.0337353,-1.89046e-05,3.91187e-09,-2.84666e-13,-69790.2,-45.8955], Tmin=(1360.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-550.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(COCsFO) + polycyclic(s2_4_4_ane) + radical(CCJCOOH)"""),
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
    label = '[O]C=C1[CH]C2(F)OOC12(1443)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,D}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-15.4154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,316,385,515,654,689,1295,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561299,0.0519637,1.81135e-05,-7.21384e-08,3.30515e-11,-1709.11,20.3715], Tmin=(100,'K'), Tmax=(989.212,'K')), NASAPolynomial(coeffs=[24.8643,0.0103474,-4.69189e-06,1.12922e-09,-9.75506e-14,-9289.27,-110.624], Tmin=(989.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.4154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + polycyclic(s2_4_4_ane) + radical(C=COJ) + radical(C=CCJCO)"""),
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
    label = '[O]C(F)C1=CC2(F)OO[C]12(1444)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {4,S} {6,S} {8,S}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-50.8544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,361,769,965,1078,1132,1246,3247,2950,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22171,0.0668096,-5.61466e-05,2.21691e-08,-3.58376e-12,-6021.87,22.2174], Tmin=(100,'K'), Tmax=(1403.57,'K')), NASAPolynomial(coeffs=[14.1562,0.0299481,-1.67527e-05,3.45789e-09,-2.50994e-13,-9652.79,-44.5692], Tmin=(1403.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.8544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CsCFHO) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(O2sj(Cs-F1sCdH)) + radical(C2CsJOO)"""),
)

species(
    label = '[O]C(F)C1(F)[CH]C2=C1OO2(1445)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {12,S}
10 C u0 p0 c0 {4,S} {8,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (135.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,2950,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.993295,0.060276,-4.6746e-05,1.65807e-08,-2.30852e-12,16379.4,29.0837], Tmin=(100,'K'), Tmax=(1686.98,'K')), NASAPolynomial(coeffs=[18.5989,0.0185311,-9.62777e-06,1.91214e-09,-1.34715e-13,10439.3,-65.0598], Tmin=(1686.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + polycyclic(s2_4_4_ene_m) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCO)"""),
)

species(
    label = '[O]C(F)=C1[CH]C2(F)OOC12(1446)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {5,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-200.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,316,385,515,654,689,1295,326,540,652,719,1357,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280427,0.0664531,-3.48302e-05,-6.88407e-09,7.22401e-12,-23982.2,22.801], Tmin=(100,'K'), Tmax=(1111.8,'K')), NASAPolynomial(coeffs=[21.6322,0.0190008,-1.04292e-05,2.24144e-09,-1.69941e-13,-30545,-90.6343], Tmin=(1111.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(CdCFO) + polycyclic(s2_4_4_ane) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = 'O[C](F)C1(F)[CH]C2(F)OOC21(1447)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {11,S} {14,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {13,S}
11 C u1 p0 c0 {3,S} {6,S} {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-398.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215464,0.093607,-0.000109691,6.84518e-08,-1.7909e-11,-47854.8,26.7783], Tmin=(100,'K'), Tmax=(905.843,'K')), NASAPolynomial(coeffs=[12.0324,0.0414248,-2.32797e-05,4.85531e-09,-3.56966e-13,-49995.6,-29.063], Tmin=(905.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-398.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)C1(F)CC2(F)OO[C]21(1448)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {9,S} {11,S}
9  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {14,S}
11 C u1 p0 c0 {5,S} {7,S} {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-380.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.768675,0.0816881,-7.97756e-05,4.15482e-08,-9.39102e-12,-45682.5,23.3785], Tmin=(100,'K'), Tmax=(1011.24,'K')), NASAPolynomial(coeffs=[10.2006,0.0443794,-2.44339e-05,5.06343e-09,-3.71087e-13,-47590,-22.2302], Tmin=(1011.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJOO)"""),
)

species(
    label = '[O][C](F)C1(F)CC2(F)OOC21(1449)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
11 C u1 p0 c0 {3,S} {6,S} {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-372.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.460423,0.0874717,-9.42776e-05,5.4062e-08,-1.31117e-11,-44716.7,24.4909], Tmin=(100,'K'), Tmax=(966.481,'K')), NASAPolynomial(coeffs=[11.5994,0.0413702,-2.2727e-05,4.7072e-09,-3.45009e-13,-46869.8,-28.8686], Tmin=(966.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-372.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(F)C1(F)[CH]C2(F)OO[C]21(1450)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {7,S} {12,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 C u1 p0 c0 {7,S} {8,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-406.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503375,0.0880072,-9.55938e-05,5.61845e-08,-1.41967e-11,-48819.5,25.7428], Tmin=(100,'K'), Tmax=(922.004,'K')), NASAPolynomial(coeffs=[10.4829,0.0447118,-2.51563e-05,5.25332e-09,-3.86619e-13,-50659.8,-21.5921], Tmin=(922.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CCJCOOH)"""),
)

species(
    label = 'FOC(F)[C]1[CH]C2(F)OOC12(1451)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {4,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
10 C u1 p0 c0 {7,S} {9,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-47.0527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,316,385,515,654,689,1295,487,638,688,1119,1325,1387,3149,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117317,0.0753493,-5.71498e-05,1.93108e-08,-2.5309e-12,-5511.13,32.8899], Tmin=(100,'K'), Tmax=(1802.7,'K')), NASAPolynomial(coeffs=[25.2573,0.0195666,-1.07339e-05,2.1455e-09,-1.50408e-13,-14575.1,-103.211], Tmin=(1802.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.0527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2sCF) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[O]C(F)[C]1C(F)C2(F)OOC12(1452)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-342.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2950,1000,259,529,569,1128,1321,1390,3140,391,562,707,872,1109,1210,1289,3137,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.235381,0.0711571,-5.04374e-05,1.50347e-08,-1.47742e-12,-41064.7,31.4019], Tmin=(100,'K'), Tmax=(1446.89,'K')), NASAPolynomial(coeffs=[21.4232,0.0232529,-1.18372e-05,2.34637e-09,-1.65766e-13,-48313,-82.5042], Tmin=(1446.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-342.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(C2CJCOOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FOC(F)C1(F)[CH][C]2OOC21(1453)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u1 p0 c0 {7,S} {10,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-28.5227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,289,311,382,485,703,1397,2750,3150,900,1100,261,493,600,1152,1365,1422,3097,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60093,0.0854099,-8.9551e-05,5.09256e-08,-1.25079e-11,-3316.97,27.1646], Tmin=(100,'K'), Tmax=(943.599,'K')), NASAPolynomial(coeffs=[10.2177,0.0446434,-2.47458e-05,5.13964e-09,-3.77211e-13,-5131.84,-18.6725], Tmin=(943.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.5227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2sCF) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CCJCOOH)"""),
)

species(
    label = '[O]C(F)C1(F)C(F)[C]2OOC21(1454)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {14,S}
11 C u1 p0 c0 {5,S} {8,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-324.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([321,327,383,445,494,1500,2950,1000,259,529,569,1128,1321,1390,3140,391,562,707,872,1109,1210,1289,3137,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751788,0.0808586,-8.16656e-05,4.51691e-08,-1.0824e-11,-38872.1,25.5565], Tmin=(100,'K'), Tmax=(964.014,'K')), NASAPolynomial(coeffs=[9.79936,0.0433177,-2.32528e-05,4.77404e-09,-3.4833e-13,-40616.5,-17.7613], Tmin=(964.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-324.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJOO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FO[CH]C1(F)[CH]C2(F)OOC12(1455)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {3,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {13,S}
11 C u1 p0 c0 {6,S} {7,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-56.8988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,289,311,382,485,703,1397,2750,3150,900,1100,316,385,515,654,689,1295,3025,407.5,1350,352.5,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0124147,0.0951738,-0.000102947,5.53367e-08,-1.21292e-11,-6705.91,26.704], Tmin=(100,'K'), Tmax=(1082.21,'K')), NASAPolynomial(coeffs=[16.2292,0.0352344,-1.98671e-05,4.15779e-09,-3.06383e-13,-10215.9,-52.8136], Tmin=(1082.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.8988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2sCF) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(CCsJO)"""),
)

species(
    label = '[O][CH]C1(F)C(F)C2(F)OOC12(1456)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {7,S} {9,S} {13,S}
11 C u1 p0 c0 {6,S} {7,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-331.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,2950,1000,333,384,448,608,1254,1480,250,417,511,1155,1315,1456,3119,3025,407.5,1350,352.5,609.728,609.728,609.729,609.729,609.732,609.732,609.732,609.733,609.734,609.736],'cm^-1')),
        HinderedRotor(inertia=(0.000453445,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0412356,0.0973378,-0.000129699,9.93394e-08,-3.22324e-11,-39764.7,26.1636], Tmin=(100,'K'), Tmax=(738.551,'K')), NASAPolynomial(coeffs=[9.51788,0.0460128,-2.54587e-05,5.24582e-09,-3.82002e-13,-41164.5,-16.6838], Tmin=(738.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-331.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + group(Cs-CsOsHH) + polycyclic(s2_4_4_ane) + radical(CCOJ) + radical(CCsJOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-150.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (148.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (155.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (14.6583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (6.46135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (6.46135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (254.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-142.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-87.4569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-103.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-52.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-150.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (109.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-47.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (149.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (14.3038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-63.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-83.4707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (220.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (220.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (129.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (237.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (58.1997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-25.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-2.11605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (6.87694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-43.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (269.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (51.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (257.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (59.7271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (237.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (94.0748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['CHFO(47)', 'FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['[O]C(F)C12[CH]C(F)(OO1)C2F(1431)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)C1(F)[CH]C2(OO2)C1F(1432)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)[C](F)C1C2OOC21F(1433)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['[O]C(F)C1(F)C2OO[C](F)C21(1434)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', 'F[CH]C1(F)[CH]C2(F)OOC12(1107)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['FC1OC2C3(F)OOC3C12F(1028)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['O=C(F)C1(F)CC2(F)OOC21(1435)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)C1(F)[CH]OOC(F)=C1(1436)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra_pri;radadd_intra_csHNd] for rate rule [Rn0c6_gamma;doublebond_intra_pri;radadd_intra_csHO]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)C(F)=CC1(F)[CH]OO1(1437)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.23261e+08,'s^-1'), n=1.06333, Ea=(126.217,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_csHNd] for rate rule [R4_Cs_RR_D;doublebond_intra_pri;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC1C(F)=CC1(F)C([O])F(1438)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.83633e+10,'s^-1'), n=0.345735, Ea=(68.5957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 68.2 to 68.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]C(F)C1=CC2(F)OOC12(1439)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(58.3508,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH]F(388)', 'FC1=CC2(F)OOC12(939)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.4589e-08,'m^3/(mol*s)'), n=3.53001, Ea=(5.06461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4944253016374622, var=1.7828810760479818, Tref=1000.0, N=135, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C(F)C1(F)C=C2OOC21(1440)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(44.0787,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=CC1(F)[CH]C2(F)OOC12(1441)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(21.6517,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHFO(47)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.66944e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(91.132,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'O=C(F)C1(F)[CH]C2(F)OOC21(1442)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(38.9568,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][CH]F(388)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F2(78)', '[O]C=C1[CH]C2(F)OOC12(1443)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(28.8413,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C(F)C1=CC2(F)OO[C]12(1444)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(245.626,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[O]C(F)C1(F)[CH]C2=C1OO2(1445)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(167.19,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]C(F)=C1[CH]C2(F)OOC12(1446)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(323.863,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['O[C](F)C1(F)[CH]C2(F)OOC21(1447)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['[O]C(F)C1(F)CC2(F)OO[C]21(1448)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.05815e+09,'s^-1'), n=0.95, Ea=(148.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['[O][C](F)C1(F)CC2(F)OOC21(1449)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.2544e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_12cy4;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['OC(F)C1(F)[CH]C2(F)OO[C]21(1450)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.76e+08,'s^-1'), n=1.2, Ea=(107.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['FOC(F)[C]1[CH]C2(F)OOC12(1451)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(100.049,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(F)[C]1C(F)C2(F)OOC12(1452)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(177.864,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['FOC(F)C1(F)[CH][C]2OOC21(1453)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(70.1813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(F)C1(F)C(F)[C]2OOC21(1454)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(167.765,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction32',
    reactants = ['FO[CH]C1(F)[CH]C2(F)OOC12(1455)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(77.9221,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O][CH]C1(F)C(F)C2(F)OOC12(1456)'],
    products = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(209.734,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #225',
    isomers = [
        '[O]C(F)C1(F)[CH]C2(F)OOC21(1023)',
    ],
    reactants = [
        ('CHFO(47)', 'FC1=CC2(F)OOC12(939)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #225',
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

