species(
    label = 'O=[C]CC(F)(F)[CH]F(7862)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {5,S} {11,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-481.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.240732,'amu*angstrom^2'), symmetry=1, barrier=(5.53491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240976,'amu*angstrom^2'), symmetry=1, barrier=(5.54052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240846,'amu*angstrom^2'), symmetry=1, barrier=(5.53752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512658,0.0889175,-0.000164636,1.54397e-07,-5.46243e-11,-57748,27.166], Tmin=(100,'K'), Tmax=(860.922,'K')), NASAPolynomial(coeffs=[7.42992,0.0300881,-1.56325e-05,3.00439e-09,-2.04872e-13,-57949.9,0.574726], Tmin=(860.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-481.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(Cds-OdCsH) + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(CCCJ=O)"""),
)

species(
    label = 'CHFCF2(55)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-511.455,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(82.003,'amu')),
        NonlinearRotor(inertia=([47.283,130.848,178.131],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([226.38,309.801,488.171,585.738,628.102,776.692,950.625,1195.27,1295.73,1386.74,1841.68,3239.84],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93201,0.00413887,7.03233e-05,-1.49827e-07,9.42397e-11,-61509.2,9.50863], Tmin=(10,'K'), Tmax=(540.182,'K')), NASAPolynomial(coeffs=[3.82032,0.0197002,-1.38027e-05,4.49179e-09,-5.49698e-13,-61712.1,7.9889], Tmin=(540.182,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-511.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=[C]CC(F)[C](F)F(7863)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-468.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.276161,'amu*angstrom^2'), symmetry=1, barrier=(6.34948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274614,'amu*angstrom^2'), symmetry=1, barrier=(6.31393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275592,'amu*angstrom^2'), symmetry=1, barrier=(6.33641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.47,'J/mol'), sigma=(5.394,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.40 K, Pc=47.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.713632,0.0822936,-0.000144918,1.33394e-07,-4.70241e-11,-56239,27.6427], Tmin=(100,'K'), Tmax=(848.531,'K')), NASAPolynomial(coeffs=[7.62081,0.0292651,-1.49936e-05,2.88885e-09,-1.98261e-13,-56674.4,-0.20403], Tmin=(848.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C(F)(F)[CH]F(3924)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {5,S} {9,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-467.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0901901,'amu*angstrom^2'), symmetry=1, barrier=(2.07365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439285,'amu*angstrom^2'), symmetry=1, barrier=(10.1,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.911684,0.0744229,-0.000105473,7.6446e-08,-2.04647e-11,-56157.1,25.9505], Tmin=(100,'K'), Tmax=(651.858,'K')), NASAPolynomial(coeffs=[10.4194,0.0250452,-1.24779e-05,2.43458e-09,-1.70888e-13,-57587.1,-17.3109], Tmin=(651.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-467.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
)

species(
    label = '[CH]F(137)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.278,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93332,-0.000263306,8.89168e-06,-1.0303e-08,3.508e-12,25853.7,4.33731], Tmin=(100,'K'), Tmax=(1056.13,'K')), NASAPolynomial(coeffs=[4.72429,0.00164127,-7.73092e-07,1.90982e-10,-1.59921e-14,25413.4,-0.815661], Tmin=(1056.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'O=[C]C[C](F)F(541)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-243.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1855,455,950,1192.27],'cm^-1')),
        HinderedRotor(inertia=(0.776172,'amu*angstrom^2'), symmetry=1, barrier=(17.8457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.61968,'amu*angstrom^2'), symmetry=1, barrier=(60.2316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83288,0.0293857,-2.14755e-05,7.65592e-09,-1.07841e-12,-29321.8,12.1615], Tmin=(10,'K'), Tmax=(1615.06,'K')), NASAPolynomial(coeffs=[9.89848,0.0143632,-7.52343e-06,1.89685e-09,-1.86956e-13,-31281.1,-20.0092], Tmin=(1615.06,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-243.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""OD[C]C[C](F)F""", comment="""Thermo library: 2-BTP_G4"""),
)

species(
    label = '[C]=O(514)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.49],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.0892,0.0020037,-1.61643e-05,2.55031e-08,-1.16411e-11,52802.7,4.52493], Tmin=(100,'K'), Tmax=(856.126,'K')), NASAPolynomial(coeffs=[0.961549,0.00569059,-3.48052e-06,7.19222e-10,-5.08058e-14,53738.7,21.4667], Tmin=(856.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C(F)(F)[CH]F(172)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 C u1 p0 c0 {4,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-330.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,1541.11],'cm^-1')),
        HinderedRotor(inertia=(0.282873,'amu*angstrom^2'), symmetry=1, barrier=(6.50381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284917,'amu*angstrom^2'), symmetry=1, barrier=(6.55079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69297,0.0567514,-9.31372e-05,8.46081e-08,-3.02383e-11,-39627.2,21.0627], Tmin=(100,'K'), Tmax=(810.463,'K')), NASAPolynomial(coeffs=[6.60547,0.0225533,-1.14234e-05,2.2402e-09,-1.56743e-13,-40096.6,0.411593], Tmin=(810.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-330.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'O=C1CC(F)(F)C1F(7973)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {4,D} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-716.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67239,0.0485781,-3.6145e-05,1.31753e-08,-1.9455e-12,-86112.7,21.0137], Tmin=(100,'K'), Tmax=(1561.12,'K')), NASAPolynomial(coeffs=[12.484,0.0208762,-9.52792e-06,1.80875e-09,-1.25262e-13,-89488.4,-35.9617], Tmin=(1561.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(Cs-(Cds-O2d)CsHH) + group(CsCCFH) + group(Cds-OdCsCs) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-CO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cyclobutanone)"""),
)

species(
    label = 'O=C=CC(F)(F)CF(7974)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-693.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24548,0.06561,-8.16372e-05,5.51099e-08,-1.5263e-11,-83367,22.3876], Tmin=(100,'K'), Tmax=(870.564,'K')), NASAPolynomial(coeffs=[9.90357,0.0258274,-1.30892e-05,2.61533e-09,-1.87704e-13,-84874.5,-18.1823], Tmin=(870.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-693.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
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
    label = 'O=[C]CC(F)=CF(6292)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {10,S}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-308.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,194,682,905,1196,1383,3221,1855,455,950,257.874],'cm^-1')),
        HinderedRotor(inertia=(0.00253413,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242922,'amu*angstrom^2'), symmetry=1, barrier=(12.1867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03987,0.0460021,-3.85314e-05,1.54489e-08,-2.53496e-12,-37042,20.3834], Tmin=(100,'K'), Tmax=(1391.97,'K')), NASAPolynomial(coeffs=[10.9506,0.0203961,-1.09382e-05,2.23347e-09,-1.61461e-13,-39522.7,-25.5527], Tmin=(1391.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=O(4373)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'F[CH][C](F)F(3610)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-285.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1767.79],'cm^-1')),
        HinderedRotor(inertia=(0.511633,'amu*angstrom^2'), symmetry=1, barrier=(11.7635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54625,0.0360263,-6.18605e-05,5.68976e-08,-2.03344e-11,-34259.7,17.1201], Tmin=(100,'K'), Tmax=(816.909,'K')), NASAPolynomial(coeffs=[5.79758,0.0130515,-6.72073e-06,1.3276e-09,-9.30772e-14,-34555.5,3.53259], Tmin=(816.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'O=C=CC(F)(F)[CH]F(7975)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u1 p0 c0 {3,S} {5,S} {10,S}
8  C u0 p0 c0 {4,D} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-492.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,3010,987.5,1337.5,450,1655,334,575,1197,1424,3202,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.230314,'amu*angstrom^2'), symmetry=1, barrier=(5.29537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.669137,'amu*angstrom^2'), symmetry=1, barrier=(15.3848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29619,0.0644542,-8.73448e-05,6.33661e-08,-1.86412e-11,-59181.4,25.2349], Tmin=(100,'K'), Tmax=(825.647,'K')), NASAPolynomial(coeffs=[9.88166,0.0228618,-1.17843e-05,2.35717e-09,-1.68776e-13,-60599.2,-14.5405], Tmin=(825.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-492.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
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
    label = 'O=[C]C=C(F)[CH]F(7976)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {4,D} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {4,S} {9,S}
7 C u1 p0 c0 {3,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-186.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.947796,'amu*angstrom^2'), symmetry=1, barrier=(21.7917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94853,'amu*angstrom^2'), symmetry=1, barrier=(21.8086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09361,0.0705841,-0.000119276,1.0604e-07,-3.73681e-11,-22286.9,21.69], Tmin=(100,'K'), Tmax=(772.986,'K')), NASAPolynomial(coeffs=[9.064,0.0224753,-1.25992e-05,2.54853e-09,-1.81273e-13,-23314.1,-13.3837], Tmin=(772.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Csj(Cd-CdF1s)(F1s)(H)) + radical(C=CCJ=O)"""),
)

species(
    label = 'CHF(40)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (138.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = '[O]C=CC(F)(F)[CH]F(7977)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u1 p0 c0 {3,S} {5,S} {10,S}
8  C u0 p0 c0 {4,S} {6,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-479.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2995,3025,975,1000,1300,1375,400,500,1630,1680,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.939138,'amu*angstrom^2'), symmetry=1, barrier=(21.5926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939947,'amu*angstrom^2'), symmetry=1, barrier=(21.6112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.606547,0.067477,-7.13939e-05,3.6236e-08,-7.10077e-12,-57514.7,27.7398], Tmin=(100,'K'), Tmax=(1254.42,'K')), NASAPolynomial(coeffs=[18.0236,0.0119381,-4.98159e-06,9.40631e-10,-6.65002e-14,-61884.4,-60.2355], Tmin=(1254.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
)

species(
    label = 'O=[C][CH]C(F)(F)CF(7978)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
7  C u1 p0 c0 {5,S} {8,S} {11,S}
8  C u1 p0 c0 {4,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-512.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.804004,0.0749453,-0.00010867,8.37572e-08,-2.57553e-11,-61565.1,24.7544], Tmin=(100,'K'), Tmax=(796.822,'K')), NASAPolynomial(coeffs=[11.1328,0.0231023,-1.10901e-05,2.12764e-09,-1.47909e-13,-63211.4,-22.7318], Tmin=(796.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C[C](F)C(F)F(7967)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-479.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735973,0.082258,-0.000146264,1.35825e-07,-4.8148e-11,-57551.4,28.2692], Tmin=(100,'K'), Tmax=(850.682,'K')), NASAPolynomial(coeffs=[7.2048,0.0299299,-1.53586e-05,2.9584e-09,-2.02882e-13,-57859.1,2.76682], Tmin=(850.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-479.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(CsCsCsF1s) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C(F)C[C](F)[CH]F(7979)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {4,D} {5,S}
8  C u1 p0 c0 {3,S} {6,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-484.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.682473,0.0860391,-0.000160831,1.55521e-07,-5.6725e-11,-58202.6,28.2075], Tmin=(100,'K'), Tmax=(849.95,'K')), NASAPolynomial(coeffs=[5.58214,0.0340639,-1.8073e-05,3.51954e-09,-2.42528e-13,-57991,11.5104], Tmin=(849.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-484.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    E0 = (-244.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-71.5858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-73.4916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (208.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (345.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-235.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-180.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-180.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (48.6109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-100.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-103.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-44.0494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (111.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-13.0028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (131.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-105.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-88.3169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-61.4717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-18.5337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CC(F)(F)[CH]F(7862)'],
    products = ['CHFCF2(55)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C]CC(F)[C](F)F(7863)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C(F)(F)[CH]F(3924)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]F(137)', 'O=[C]C[C](F)F(541)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=O(514)', '[CH2]C(F)(F)[CH]F(172)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]CC(F)(F)[CH]F(7862)'],
    products = ['O=C1CC(F)(F)C1F(7973)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_1H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=[C]CC(F)(F)[CH]F(7862)'],
    products = ['O=C=CC(F)(F)CF(7974)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CO(13)', '[CH2]C(F)(F)[CH]F(172)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(31.0376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'O=[C]CC(F)=CF(6292)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.2858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CHFCF2(55)', '[CH2][C]=O(4373)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(13.5501,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH][C](F)F(3610)', 'CH2CO(28)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(5.28467,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C=CC(F)(F)[CH]F(7975)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C](F)F(3610)', '[CH2][C]=O(4373)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -4.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'O=[C]C=C(F)[CH]F(7976)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3027.76,'m^3/(mol*s)'), n=0.596786, Ea=(217.251,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.06681560721952781, var=6.699388179313154, Tref=1000.0, N=4, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHF(40)', 'O=[C]C[C](F)F(541)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C=CC(F)(F)[CH]F(7977)'],
    products = ['O=[C]CC(F)(F)[CH]F(7862)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;XH_out] for rate rule [R2H_S;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]CC(F)(F)[CH]F(7862)'],
    products = ['O=[C][CH]C(F)(F)CF(7978)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(182.734,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]CC(F)(F)[CH]F(7862)'],
    products = ['O=[C]C[C](F)C(F)F(7967)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(182.638,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C]CC(F)(F)[CH]F(7862)'],
    products = ['O=C(F)C[C](F)[CH]F(7979)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(225.576,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #2086',
    isomers = [
        'O=[C]CC(F)(F)[CH]F(7862)',
    ],
    reactants = [
        ('CHFCF2(55)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2086',
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

