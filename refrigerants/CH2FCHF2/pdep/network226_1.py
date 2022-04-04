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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700842,0.0786221,-6.8804e-05,2.94181e-08,-5.20195e-12,-44161.7,25.5463], Tmin=(100,'K'), Tmax=(1294.88,'K')), NASAPolynomial(coeffs=[14.5168,0.0359432,-1.93642e-05,3.96409e-09,-2.8757e-13,-47739.7,-44.6779], Tmin=(1294.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsCsF1s)"""),
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
    label = '[O]C(F)C1[C](F)C(F)C12OO2(1397)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {1,S} {7,S} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {6,S} {8,S} {14,S}
11 C u1 p0 c0 {3,S} {8,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-360.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,259,529,569,1128,1321,1390,3140,391,562,707,872,1109,1210,1289,3137,212,367,445,1450,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744243,0.0915948,-8.65449e-05,3.42499e-08,-3.77441e-12,-43202.4,32.3112], Tmin=(100,'K'), Tmax=(1059.27,'K')), NASAPolynomial(coeffs=[24.5312,0.0156382,-6.58157e-06,1.29215e-09,-9.51527e-14,-49650.4,-96.2442], Tmin=(1059.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s1_3_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)[CH]C1(F)C2OOC21F(1398)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {7,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {11,S} {13,S}
11 C u1 p0 c0 {7,S} {10,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-357.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,333,384,448,608,1254,1480,2950,1000,391,562,707,872,1109,1210,1289,3137,3025,407.5,1350,352.5,543.603,543.746,543.915,544.035,544.039,544.096,544.109,544.327],'cm^-1')),
        HinderedRotor(inertia=(0.000571198,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000570258,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184797,0.0840247,-7.87741e-05,3.54917e-08,-6.4068e-12,-42875.8,29.2121], Tmin=(100,'K'), Tmax=(1310.58,'K')), NASAPolynomial(coeffs=[18.4272,0.0283474,-1.50497e-05,3.07643e-09,-2.23414e-13,-47657.4,-63.7309], Tmin=(1310.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + polycyclic(s2_3_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

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
    label = '[O]C(F)C1C2(F)[CH]OOC12F(1399)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {9,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {13,S}
11 C u1 p0 c0 {5,S} {8,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-420.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.500192,0.079934,-7.10173e-05,3.09992e-08,-5.53434e-12,-50497.1,26.7702], Tmin=(100,'K'), Tmax=(1301.31,'K')), NASAPolynomial(coeffs=[15.6448,0.0333823,-1.73578e-05,3.50914e-09,-2.53113e-13,-54438.7,-50.2822], Tmin=(1301.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFHO) + polycyclic(s2_3_5_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'F[CH]C1[C](F)C2OOC12F(1066)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u1 p0 c0 {3,S} {7,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-205.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,334,575,1197,1424,3202,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28633,0.0692787,-5.77154e-05,2.35263e-08,-4.07237e-12,-24578,21.3296], Tmin=(100,'K'), Tmax=(1279.26,'K')), NASAPolynomial(coeffs=[11.3603,0.03778,-2.07823e-05,4.27954e-09,-3.11153e-13,-27155.5,-29.7526], Tmin=(1279.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsFHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'FC1OC2(F)C3OOC3(F)C12(1027)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {5,S} {10,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
10 C u0 p0 c0 {6,S} {8,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-619.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00659773,0.0764496,-5.42648e-05,1.56243e-08,-1.50761e-12,-74354.7,18.1848], Tmin=(100,'K'), Tmax=(1558.45,'K')), NASAPolynomial(coeffs=[25.9367,0.0214986,-1.2575e-05,2.58164e-09,-1.8395e-13,-83854.1,-123.021], Tmin=(1558.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) - ring(Cyclobutane)"""),
)

species(
    label = 'OC(F)C1=C(F)C2OOC12F(1400)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {8,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-673.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16658,0.0865777,-8.2411e-05,3.74742e-08,-6.84959e-12,-80830.7,24.4705], Tmin=(100,'K'), Tmax=(1289.75,'K')), NASAPolynomial(coeffs=[18.4244,0.0299528,-1.65546e-05,3.43292e-09,-2.51081e-13,-85540.3,-68.2582], Tmin=(1289.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-673.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'O=C(F)C1C(F)C2OOC21F(1401)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {7,S} {9,S} {14,S}
11 C u0 p0 c0 {3,S} {6,D} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-733.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607,0.0765938,-6.25469e-05,2.39528e-08,-3.70345e-12,-88057.4,22.4979], Tmin=(100,'K'), Tmax=(1485.47,'K')), NASAPolynomial(coeffs=[17.4522,0.031234,-1.67436e-05,3.39669e-09,-2.43935e-13,-93062,-65.4367], Tmin=(1485.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(COCsFO) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'OC(F)C1C(F)=C2OOC21F(1030)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
10 C u0 p0 c0 {5,S} {8,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-602.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785208,0.0811429,-8.26055e-05,4.51679e-08,-1.06806e-11,-72297.6,23.2736], Tmin=(100,'K'), Tmax=(974.505,'K')), NASAPolynomial(coeffs=[10.0879,0.0429572,-2.38259e-05,4.95466e-09,-3.6389e-13,-74110.6,-21.3657], Tmin=(974.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-602.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(CsCFHO) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane)"""),
)

species(
    label = '[O]C(F)C1[C](F)OOC=C1F(1402)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {10,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {11,D}
10 C u1 p0 c0 {3,S} {4,S} {7,S}
11 C u0 p0 c0 {5,S} {9,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-430.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,391,562,707,872,1109,1210,1289,3137,323,467,575,827,1418,395,473,707,1436,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.180154,0.086163,-9.77564e-05,5.91138e-08,-1.44222e-11,-51605.6,29.1472], Tmin=(100,'K'), Tmax=(992.541,'K')), NASAPolynomial(coeffs=[14.0015,0.0304609,-1.35734e-05,2.5688e-09,-1.79405e-13,-54349.2,-37.429], Tmin=(992.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-430.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCFHO) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C(F)C=C(F)C1OO[C]1F(1403)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 C u1 p0 c0 {3,S} {5,S} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-379.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,361,769,965,1078,1132,1246,3247,323,467,575,827,1418,3010,987.5,1337.5,450,1655,395,473,707,1436,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.350232,0.085486,-9.69982e-05,5.75504e-08,-1.39201e-11,-45522.4,35.1554], Tmin=(100,'K'), Tmax=(991.808,'K')), NASAPolynomial(coeffs=[13.5911,0.0320841,-1.62326e-05,3.26113e-09,-2.35433e-13,-48148.9,-28.6151], Tmin=(991.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(O2sj(Cs-F1sCdH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC1(F)C=C(F)C1C([O])F(1404)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {9,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {7,S} {13,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-456.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,323,467,575,827,1418,180,180,180,180,431.223,1022.09,1022.52,1895.96],'cm^-1')),
        HinderedRotor(inertia=(0.75504,'amu*angstrom^2'), symmetry=1, barrier=(17.3599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95542,'amu*angstrom^2'), symmetry=1, barrier=(44.959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0118136,0.0899124,-0.000101018,5.78137e-08,-1.32438e-11,-54814.7,31.0992], Tmin=(100,'K'), Tmax=(1055.44,'K')), NASAPolynomial(coeffs=[16.203,0.0285491,-1.38076e-05,2.72651e-09,-1.95312e-13,-58232.4,-47.8872], Tmin=(1055.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sH)) + radical(ROOJ)"""),
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
    label = '[O]C(F)C1=C(F)C2OOC12F(1405)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {8,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-445.976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,361,769,965,1078,1132,1246,3247,323,467,575,827,1418,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282413,0.0813388,-7.54041e-05,3.27596e-08,-5.66402e-12,-53504.8,24.8972], Tmin=(100,'K'), Tmax=(1367.16,'K')), NASAPolynomial(coeffs=[19.3129,0.0256604,-1.43163e-05,2.97166e-09,-2.17024e-13,-58708.4,-72.8656], Tmin=(1367.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(O2sj(Cs-F1sCdH))"""),
)

species(
    label = '[O]C(F)C1C(F)=C2OOC21F(1406)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
10 C u0 p0 c0 {5,S} {8,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-375.433,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,323,467,575,827,1418,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82767,0.0795346,-9.10459e-05,5.93669e-08,-1.68561e-11,-45048.4,24.0298], Tmin=(100,'K'), Tmax=(825.253,'K')), NASAPolynomial(coeffs=[8.4787,0.0424503,-2.3641e-05,4.91533e-09,-3.60823e-13,-46311.2,-11.4126], Tmin=(825.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-375.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFO) + group(CsCFHO) + group(Cds-CdsCsOs) + group(CdCsCdF) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane) + radical(O2sj(Cs-CsF1sH)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    label = 'O=CC1[C](F)C2OOC12F(1407)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {6,S} {9,S} {12,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {5,D} {7,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-288.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4285,0.065725,-5.16583e-05,1.94861e-08,-3.09847e-12,-34635.8,21.3226], Tmin=(100,'K'), Tmax=(1376.21,'K')), NASAPolynomial(coeffs=[11.4988,0.0364552,-1.97558e-05,4.03181e-09,-2.91071e-13,-37407.6,-30.4767], Tmin=(1376.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s)"""),
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
    label = 'O=C(F)C1[C](F)C2OOC21F(1408)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {5,S} {7,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {8,S} {9,S}
11 C u0 p0 c0 {3,S} {6,D} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-542.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2750,3150,900,1100,212,367,445,1450,486,617,768,1157,1926,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997384,0.0761855,-7.39797e-05,3.76901e-08,-8.29347e-12,-65168.5,23.8321], Tmin=(100,'K'), Tmax=(1037.18,'K')), NASAPolynomial(coeffs=[10.1948,0.0407157,-2.26835e-05,4.71941e-09,-3.46476e-13,-67076.5,-20.8759], Tmin=(1037.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-542.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(COCsFO) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s)"""),
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
    label = '[O]C(F)C1=C(F)C2OO[C]12(1409)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {4,S} {6,S} {8,S}
10 C u0 p0 c0 {2,S} {6,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-14.7065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,361,769,965,1078,1132,1246,3247,323,467,575,827,1418,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.810749,0.0658999,-5.38413e-05,2.01294e-08,-2.96495e-12,-1651,26.9895], Tmin=(100,'K'), Tmax=(1594.43,'K')), NASAPolynomial(coeffs=[18.9367,0.0204267,-1.10613e-05,2.24216e-09,-1.60307e-13,-7431.12,-68.914], Tmin=(1594.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.7065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(O2sj(Cs-F1sCdH)) + radical(C2CsJOO)"""),
)

species(
    label = '[O]C=C1[C](F)C2OOC12F(1410)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {4,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-195.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,346,659,817,1284,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145357,0.0696706,-4.31998e-05,-1.00761e-09,6.0664e-12,-23339.2,22.9589], Tmin=(100,'K'), Tmax=(1085.8,'K')), NASAPolynomial(coeffs=[22.3701,0.0162118,-8.60316e-06,1.85276e-09,-1.41657e-13,-29840.5,-93.8052], Tmin=(1085.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + polycyclic(s2_4_4_ane) + radical(C=COJ) + radical(CsCdCsF1s)"""),
)

species(
    label = '[O]C(F)C1[C](F)C2=C1OO2(1411)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u1 p0 c0 {2,S} {6,S} {10,S}
10 C u0 p0 c0 {4,S} {8,D} {9,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (96.3979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,391,562,707,872,1109,1210,1289,3137,346,659,817,1284,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.871273,0.0574815,-3.29023e-05,8.82513e-10,3.08981e-12,11716,28.4349], Tmin=(100,'K'), Tmax=(1174.36,'K')), NASAPolynomial(coeffs=[17.3474,0.0206517,-1.04983e-05,2.15082e-09,-1.57684e-13,6516.11,-59.3641], Tmin=(1174.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.3979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + polycyclic(s2_4_4_ene_m) + radical(O2sj(Cs-CsF1sH)) + radical(CsCdCsF1s)"""),
)

species(
    label = '[O]C(F)[C]1C(F)C2OOC12F(1412)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {9,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-367.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.239074,0.0723908,-5.34592e-05,1.76764e-08,-2.26832e-12,-44078.4,30.0879], Tmin=(100,'K'), Tmax=(1839.29,'K')), NASAPolynomial(coeffs=[24.598,0.0194171,-1.02582e-05,2.01807e-09,-1.4005e-13,-53039.2,-102.274], Tmin=(1839.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(C2CJCOOH)"""),
)

species(
    label = '[O]C(F)C1C(F)[C]2OOC21F(1413)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {14,S}
11 C u1 p0 c0 {5,S} {8,S} {9,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-371.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,316,385,515,654,689,1295,259,529,569,1128,1321,1390,3140,391,562,707,872,1109,1210,1289,3137,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.841903,0.0771602,-6.57306e-05,2.74628e-08,-4.79974e-12,-44631.3,23.1344], Tmin=(100,'K'), Tmax=(1292.59,'K')), NASAPolynomial(coeffs=[13.3837,0.0383487,-2.0691e-05,4.23302e-09,-3.06838e-13,-47873.5,-40.5911], Tmin=(1292.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-371.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJOO)"""),
)

species(
    label = 'O[C](F)C1[C](F)C2OOC21F(1414)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {11,S} {14,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {7,S} {9,S}
11 C u1 p0 c0 {3,S} {6,S} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-400.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375616,0.0898897,-0.000102214,6.26559e-08,-1.62525e-11,-48004,26.6525], Tmin=(100,'K'), Tmax=(908.836,'K')), NASAPolynomial(coeffs=[11.221,0.0421578,-2.34365e-05,4.87067e-09,-3.57508e-13,-49975.4,-24.6336], Tmin=(908.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-400.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O][C](F)C1C(F)C2OOC21F(1415)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {7,S} {9,S} {14,S}
11 C u1 p0 c0 {3,S} {6,S} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-364.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,333,384,448,608,1254,1480,250,417,511,1155,1315,1456,3119,395,473,707,1436,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549763,0.0826381,-7.87752e-05,3.77703e-08,-7.51906e-12,-43665.7,24.1989], Tmin=(100,'K'), Tmax=(1166.54,'K')), NASAPolynomial(coeffs=[13.9245,0.036777,-1.98047e-05,4.06921e-09,-2.96632e-13,-46786.2,-42.3866], Tmin=(1166.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(F)[C]1[C](F)C2OOC12F(1416)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
10 C u1 p0 c0 {7,S} {9,S} {11,S}
11 C u1 p0 c0 {3,S} {8,S} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-403.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809518,0.07198,-5.45935e-05,1.87377e-08,-2.56356e-12,-48452.1,29.786], Tmin=(100,'K'), Tmax=(1663.15,'K')), NASAPolynomial(coeffs=[18.6857,0.028987,-1.58184e-05,3.195e-09,-2.27262e-13,-54398.3,-65.5503], Tmin=(1663.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(C2CJCOOH) + radical(CsCsCsF1s)"""),
)

species(
    label = 'OC(F)C1[C](F)[C]2OOC21F(1417)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {9,S} {14,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u1 p0 c0 {3,S} {7,S} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-408.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667036,0.0842537,-8.80136e-05,5.02815e-08,-1.2504e-11,-48969,25.604], Tmin=(100,'K'), Tmax=(928.772,'K')), NASAPolynomial(coeffs=[9.68357,0.0454216,-2.52986e-05,5.26504e-09,-3.86845e-13,-50643.8,-17.2293], Tmin=(928.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-408.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FOC(F)C1[C](F)C2OO[C]21(1418)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {4,S} {10,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 C u1 p0 c0 {2,S} {7,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-29.7171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,261,493,600,1152,1365,1422,3097,212,367,445,1450,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767435,0.0816291,-8.19025e-05,4.49661e-08,-1.08025e-11,-3466.51,27.0152], Tmin=(100,'K'), Tmax=(954.837,'K')), NASAPolynomial(coeffs=[9.43396,0.045324,-2.48701e-05,5.14691e-09,-3.77056e-13,-5121.57,-14.3954], Tmin=(954.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.7171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2sCF) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(F)C1[C]2OOC2C1(F)F(1419)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {4,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {6,S} {7,S} {14,S}
11 C u1 p0 c0 {5,S} {7,S} {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-360.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,222,329,445,522,589,1214,1475,391,562,707,872,1109,1210,1289,3137,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854568,0.0773906,-6.96789e-05,3.21554e-08,-6.3069e-12,-43202.1,25.1559], Tmin=(100,'K'), Tmax=(1160.35,'K')), NASAPolynomial(coeffs=[11.7155,0.0399507,-2.12803e-05,4.34874e-09,-3.15942e-13,-45722.6,-28.8572], Tmin=(1160.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + polycyclic(s2_4_4_ane) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJOO)"""),
)

species(
    label = 'FO[CH]C1[C](F)C2OOC12F(1420)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {3,S} {11,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {7,S} {9,S}
11 C u1 p0 c0 {6,S} {7,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-50.7569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,333,384,448,608,1254,1480,212,367,445,1450,3025,407.5,1350,352.5,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.154973,0.0906217,-9.25806e-05,4.67129e-08,-9.61216e-12,-5971.28,26.8194], Tmin=(100,'K'), Tmax=(1146.69,'K')), NASAPolynomial(coeffs=[16.2725,0.0343991,-1.90355e-05,3.95513e-09,-2.9021e-13,-9667.65,-53.1445], Tmin=(1146.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.7569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2sCF) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsOsHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(CCsJO)"""),
)

species(
    label = '[O][CH]C1C(F)(F)C2OOC12F(1421)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
11 C u1 p0 c0 {6,S} {7,S} {14,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-385.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333739,0.0914149,-0.000105856,6.78494e-08,-1.86161e-11,-46232.4,23.9918], Tmin=(100,'K'), Tmax=(858.468,'K')), NASAPolynomial(coeffs=[10.2348,0.0452805,-2.52449e-05,5.24747e-09,-3.85085e-13,-47932.3,-22.2641], Tmin=(858.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsOsHH) + polycyclic(s2_4_4_ane) + radical(CCOJ) + radical(CCsJOH)"""),
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
    E0 = (-122.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (184.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (45.5397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (36.2394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (37.6618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (283.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-113.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-58.8732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-58.8732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-97.3001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-53.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-8.53846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-122.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (3.24923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (26.5985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (90.6513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (49.5049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-34.1113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-49.9772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (250.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (178.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (89.1283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (242.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (51.7813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-2.21931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (3.18672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (13.8184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-46.9927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-89.5363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (286.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (109.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (271.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (97.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['CHFO(47)', 'FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(F)C1[C](F)C(F)C12OO2(1397)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)[CH]C1(F)C2OOC21F(1398)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)C1(F)[CH]C2(F)OOC21(1023)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['[O]C(F)C1C2(F)[CH]OOC12F(1399)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', 'F[CH]C1[C](F)C2OOC12F(1066)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['FC1OC2(F)C3OOC3(F)C12(1027)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['OC(F)C1=C(F)C2OOC12F(1400)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['O=C(F)C1C(F)C2OOC21F(1401)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['OC(F)C1C(F)=C2OOC21F(1030)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)C1[C](F)OOC=C1F(1402)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c6_gamma;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)C=C(F)C1OO[C]1F(1403)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC1(F)C=C(F)C1C([O])F(1404)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.10864e+11,'s^-1'), n=0.310877, Ea=(88.8126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 88.0 to 88.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH]F(388)', 'FC1=CC2(F)OOC12(939)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.43165,'m^3/(mol*s)'), n=1.58832, Ea=(25.6817,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.2743586997156029, var=2.1590031255329776, Tref=1000.0, N=359, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[O]C(F)C1=C(F)C2OOC12F(1405)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(205.9,'m^3/(mol*s)'), n=1.596, Ea=(14.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]C(F)C1C(F)=C2OOC21F(1406)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(45.45,'m^3/(mol*s)'), n=1.71, Ea=(8.43039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_N-Sp-9R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_N-Sp-9R!H=5R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', 'O=CC1[C](F)C2OOC12F(1407)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(19.4421,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHFO(47)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.66944e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(90.2567,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'O=C(F)C1[C](F)C2OOC21F(1408)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(35.0398,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][CH]F(388)', 'F[C]1[CH]C2(F)OOC12(943)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C(F)C1=C(F)C2OO[C]12(1409)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(228.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[O]C=C1[C](F)C2OOC12F(1410)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(319.707,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]C(F)C1[C](F)C2=C1OO2(1411)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(181.713,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['[O]C(F)[C]1C(F)C2OOC12F(1412)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.93363e+09,'s^-1'), n=1.033, Ea=(174.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_Cs2] for rate rule [R2H_S_cy4;C_rad_out_noH;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)C1C(F)[C]2OOC21F(1413)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.26888e+09,'s^-1'), n=0.972714, Ea=(123.905,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_NDMustO;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['O[C](F)C1[C](F)C2OOC21F(1414)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C](F)C1C(F)C2OOC21F(1415)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(17992.9,'s^-1'), n=2.5106, Ea=(132.012,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_23cy4;Y_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['OC(F)[C]1[C](F)C2OOC12F(1416)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['OC(F)C1[C](F)[C]2OOC21F(1417)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(214655,'s^-1'), n=1.70206, Ea=(32.737,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_CCC;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['FOC(F)C1[C](F)C2OO[C]21(1418)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(70.1813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(F)C1[C]2OOC2C1(F)F(1419)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(223.716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['FO[CH]C1[C](F)C2OOC12F(1420)'],
    products = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(75.9345,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C(F)C1[C](F)C2OOC21F(1024)'],
    products = ['[O][CH]C1C(F)(F)C2OOC12F(1421)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(219.941,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #226',
    isomers = [
        '[O]C(F)C1[C](F)C2OOC21F(1024)',
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
    label = 'PDepNetwork #226',
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

