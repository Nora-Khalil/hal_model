species(
    label = 'F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {8,S} {9,S} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1145.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,350,440,435,1725,3025,407.5,1350,352.5,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,1952.8],'cm^-1')),
        HinderedRotor(inertia=(0.260763,'amu*angstrom^2'), symmetry=1, barrier=(5.99545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260352,'amu*angstrom^2'), symmetry=1, barrier=(5.986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262693,'amu*angstrom^2'), symmetry=1, barrier=(6.03984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75801,'amu*angstrom^2'), symmetry=1, barrier=(40.42,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50681,0.133506,-0.000219857,1.94354e-07,-6.84139e-11,-137560,41.5557], Tmin=(100,'K'), Tmax=(775.942,'K')), NASAPolynomial(coeffs=[12.9424,0.0458445,-2.49268e-05,4.99304e-09,-3.53626e-13,-139406,-21.9322], Tmin=(775.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1145.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Cs_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'FC=C=C(F)F(1325)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-364.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,94,120,354,641,825,1294,540,610,2055,2775.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.43,'J/mol'), sigma=(4.61545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.73 K, Pc=67.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78631,0.0191407,1.84311e-05,-5.99612e-08,3.72113e-11,-43861.2,10.8341], Tmin=(10,'K'), Tmax=(620.754,'K')), NASAPolynomial(coeffs=[5.19832,0.0213405,-1.41861e-05,4.38954e-09,-5.1374e-13,-44254.2,2.94178], Tmin=(620.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-364.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FCDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=CC(F)F(344)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u0 p0 c0 {5,S} {7,D} {9,S}
7 C u0 p0 c0 {3,S} {4,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-792.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.142102,'amu*angstrom^2'), symmetry=1, barrier=(3.2672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2858.21,'J/mol'), sigma=(4.57471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=446.45 K, Pc=67.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74532,0.0233589,0.000133412,-5.9431e-07,7.11464e-10,-95378,11.8568], Tmin=(10,'K'), Tmax=(302.26,'K')), NASAPolynomial(coeffs=[4.97982,0.0307765,-2.12832e-05,6.89182e-09,-8.42907e-13,-95561.1,5.58312], Tmin=(302.26,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-792.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FC(F)DCC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(=C(F)F)C([C](F)F)C(F)F(6742)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {15,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {3,S} {4,S} {8,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1124.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,190,488,555,1236,1407,234,589,736,816,1240,3237,182,240,577,636,1210,1413,344,345.27,2452.52],'cm^-1')),
        HinderedRotor(inertia=(0.119488,'amu*angstrom^2'), symmetry=1, barrier=(10.0644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118097,'amu*angstrom^2'), symmetry=1, barrier=(10.061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447531,'amu*angstrom^2'), symmetry=1, barrier=(37.8864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.451329,'amu*angstrom^2'), symmetry=1, barrier=(37.9073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3373.95,'J/mol'), sigma=(5.28602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.00 K, Pc=51.83 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.40762,0.128691,-0.000190344,1.47563e-07,-4.57753e-11,-135004,41.3091], Tmin=(100,'K'), Tmax=(788.372,'K')), NASAPolynomial(coeffs=[16.1634,0.0395426,-2.07292e-05,4.13617e-09,-2.94354e-13,-137774,-39.283], Tmin=(788.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1124.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsFFH) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(CdCFF) + radical(CsCsF1sF1s) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {5,S} {11,S} {15,S}
11 C u1 p0 c0 {8,S} {10,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {13,D}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1021.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,180,180,2241.89],'cm^-1')),
        HinderedRotor(inertia=(0.203212,'amu*angstrom^2'), symmetry=1, barrier=(4.67223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203187,'amu*angstrom^2'), symmetry=1, barrier=(4.67166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204125,'amu*angstrom^2'), symmetry=1, barrier=(4.69323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55209,'amu*angstrom^2'), symmetry=1, barrier=(35.6856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3348.87,'J/mol'), sigma=(5.57238,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=523.09 K, Pc=43.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32681,0.129436,-0.000211359,1.88309e-07,-6.72684e-11,-122624,42.6331], Tmin=(100,'K'), Tmax=(763.254,'K')), NASAPolynomial(coeffs=[11.9671,0.0474391,-2.59866e-05,5.23497e-09,-3.72467e-13,-124294,-15.5587], Tmin=(763.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1021.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsCsHH) + group(CsCCFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'F[CH]C(=C(F)[CH]C(F)F)C(F)(F)F(7387)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {9,S} {12,S} {15,S}
12 C u0 p0 c0 {6,S} {10,D} {11,S}
13 C u1 p0 c0 {7,S} {10,S} {16,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1235.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.747139,0.111708,-0.000136926,8.72277e-08,-2.24775e-11,-148437,36.4618], Tmin=(100,'K'), Tmax=(937.746,'K')), NASAPolynomial(coeffs=[16.4086,0.038527,-1.98622e-05,4.00139e-09,-2.88842e-13,-151655,-45.2015], Tmin=(937.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1235.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cds(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = '[CH]C(F)F-2(957)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u2 p0 c0 {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-68.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,200.056,1350.78,1616.29],'cm^-1')),
        HinderedRotor(inertia=(0.00421206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96601,0.0195049,-1.18279e-05,2.09006e-09,3.19831e-13,-8190.14,12.8831], Tmin=(100,'K'), Tmax=(1211.04,'K')), NASAPolynomial(coeffs=[7.87163,0.00800485,-3.40879e-06,6.6204e-10,-4.73288e-14,-9723.19,-13.1468], Tmin=(1211.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = 'F[CH]C([C](F)F)=C(F)F(3166)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {7,S} {8,S} {9,D}
7  C u1 p0 c0 {1,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-683.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,234,589,736,816,1240,3237,161,297,490,584,780,1358,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.000412985,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000413036,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29869,0.0651463,-8.30537e-05,5.7336e-08,-1.62574e-11,-82156.9,26.972], Tmin=(100,'K'), Tmax=(850.021,'K')), NASAPolynomial(coeffs=[9.68445,0.0256869,-1.34247e-05,2.72914e-09,-1.97813e-13,-83582.5,-12.1222], Tmin=(850.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CsCd)(F1s)(H)) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = '[CH]F(804)',
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
    label = 'FC(F)=[C]C(F)(F)[CH]C(F)F(7388)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {6,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-808.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,562,600,623,1070,1265,1685,370,180,180,1916.04],'cm^-1')),
        HinderedRotor(inertia=(0.384639,'amu*angstrom^2'), symmetry=1, barrier=(8.8436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385501,'amu*angstrom^2'), symmetry=1, barrier=(8.86342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384492,'amu*angstrom^2'), symmetry=1, barrier=(8.84022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.458849,0.112452,-0.0002039,1.93158e-07,-7.05649e-11,-97134.5,36.5504], Tmin=(100,'K'), Tmax=(816.964,'K')), NASAPolynomial(coeffs=[8.47346,0.0413517,-2.31113e-05,4.62783e-09,-3.25547e-13,-97680.7,0.852136], Tmin=(816.964,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-808.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Cs_S) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sF1s))"""),
)

species(
    label = 'FC(F)=C1C(F)C(C(F)F)C1(F)F(7271)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {8,S} {12,S} {15,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {8,S} {16,S}
12 C u0 p0 c0 {9,S} {10,S} {13,D}
13 C u0 p0 c0 {6,S} {7,S} {12,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1347.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.489668,0.109507,-0.000140147,9.90058e-08,-2.91537e-11,-161937,30.0363], Tmin=(100,'K'), Tmax=(814.784,'K')), NASAPolynomial(coeffs=[12.1918,0.0472497,-2.55333e-05,5.22788e-09,-3.79896e-13,-164003,-28.547], Tmin=(814.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1347.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FCC(=C(F)F)C(F)(F)C=C(F)F(6751)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u0 p0 c0 {8,S} {13,D} {16,S}
12 C u0 p0 c0 {4,S} {5,S} {10,D}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1341.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.40836,0.130916,-0.000214107,1.88278e-07,-6.58461e-11,-161154,37.4374], Tmin=(100,'K'), Tmax=(783.004,'K')), NASAPolynomial(coeffs=[12.8055,0.0450092,-2.40676e-05,4.78956e-09,-3.3793e-13,-162973,-25.0568], Tmin=(783.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1341.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF)"""),
)

species(
    label = 'FC(F)[CH]C(F)(F)[C]1C(F)C1(F)F(7389)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {9,S} {12,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
11 C u0 p0 c0 {6,S} {7,S} {13,S} {15,S}
12 C u1 p0 c0 {8,S} {9,S} {10,S}
13 C u1 p0 c0 {10,S} {11,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1057.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24311,0.12853,-0.000213622,1.92881e-07,-6.86972e-11,-126997,40.6957], Tmin=(100,'K'), Tmax=(803.197,'K')), NASAPolynomial(coeffs=[10.7521,0.0483147,-2.55726e-05,5.05383e-09,-3.54575e-13,-128263,-10.4332], Tmin=(803.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1057.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsFFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH][C]1C(F)(F)C(C(F)F)C1(F)F(7390)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {8,S} {15,S}
12 C u1 p0 c0 {9,S} {10,S} {13,S}
13 C u1 p0 c0 {7,S} {12,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1156.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194032,0.109322,-0.000119314,1.22018e-08,5.43098e-11,-139014,33.9459], Tmin=(100,'K'), Tmax=(502.144,'K')), NASAPolynomial(coeffs=[11.1189,0.0486572,-2.60779e-05,5.22556e-09,-3.71387e-13,-140521,-16.5369], Tmin=(502.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1156.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsFFH) + group(CsCsFHH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1([C](F)F)C(C(F)F)C1(F)F(7317)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {8,S} {10,S} {11,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {9,S} {15,S}
12 C u1 p0 c0 {5,S} {8,S} {16,S}
13 C u1 p0 c0 {6,S} {7,S} {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1053.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10606,0.120153,-0.000152376,9.81859e-08,-2.53609e-11,-126564,36.7687], Tmin=(100,'K'), Tmax=(938.998,'K')), NASAPolynomial(coeffs=[18.3822,0.0371355,-1.97601e-05,4.03115e-09,-2.92913e-13,-130224,-56.024], Tmin=(938.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1053.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsFHH) + group(CsCsFFH) + group(CsCsFFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[CH]C(C(F)=CC(F)F)=C(F)F(7391)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {13,S}
8  C u0 p0 c0 {10,S} {11,S} {12,D}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {8,S} {9,D}
11 C u1 p0 c0 {4,S} {8,S} {15,S}
12 C u0 p0 c0 {5,S} {6,S} {8,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-996.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,3010,987.5,1337.5,450,1655,280,518,736,852,873,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,2427.27],'cm^-1')),
        HinderedRotor(inertia=(0.502656,'amu*angstrom^2'), symmetry=1, barrier=(11.5571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41474,'amu*angstrom^2'), symmetry=1, barrier=(32.5277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41517,'amu*angstrom^2'), symmetry=1, barrier=(32.5375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.603041,0.11017,-0.000148231,1.05306e-07,-3.03198e-11,-119707,33.5982], Tmin=(100,'K'), Tmax=(842.564,'K')), NASAPolynomial(coeffs=[14.5276,0.0383387,-2.03513e-05,4.12347e-09,-2.97563e-13,-122256,-36.8063], Tmin=(842.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-996.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF) + radical(CsCdF1sH)"""),
)

species(
    label = 'F[CH][C]=C(F)F(2842)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-200.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.793],'cm^-1')),
        HinderedRotor(inertia=(0.35565,'amu*angstrom^2'), symmetry=1, barrier=(8.17708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377452,-4.40203e-05,2.68135e-08,-6.58286e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.929,'K')), NASAPolynomial(coeffs=[8.46193,0.0130101,-6.31222e-06,1.26456e-09,-9.14162e-14,-25275.2,-12.1012], Tmin=(983.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)(F)C=CF(7392)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {7,S} {12,D} {13,S}
10 C u1 p0 c0 {4,S} {8,S} {15,S}
11 C u0 p0 c0 {5,S} {6,S} {8,D}
12 C u0 p0 c0 {3,S} {9,D} {14,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-994.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.312,'amu*angstrom^2'), symmetry=1, barrier=(30.1654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31284,'amu*angstrom^2'), symmetry=1, barrier=(30.1848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31389,'amu*angstrom^2'), symmetry=1, barrier=(30.209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (189.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.619489,0.112775,-0.000161706,1.18758e-07,-3.20705e-11,-119471,34.7972], Tmin=(100,'K'), Tmax=(625.999,'K')), NASAPolynomial(coeffs=[12.8283,0.0405362,-2.14122e-05,4.28296e-09,-3.04983e-13,-121423,-25.9241], Tmin=(625.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-994.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFH) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[CH]C(=C(F)F)C(F)(F)C=C(F)F(7393)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {13,S}
4  F u0 p3 c0 {13,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {11,S} {12,D}
10 C u0 p0 c0 {8,S} {13,D} {14,S}
11 C u1 p0 c0 {5,S} {9,S} {15,S}
12 C u0 p0 c0 {6,S} {7,S} {9,D}
13 C u0 p0 c0 {3,S} {4,S} {10,D}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-1195.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,141,223,164,316,539,615,578,694,1133,1287,1372,1454,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32311,'amu*angstrom^2'), symmetry=1, barrier=(30.4208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31776,'amu*angstrom^2'), symmetry=1, barrier=(30.298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31614,'amu*angstrom^2'), symmetry=1, barrier=(30.2607,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (207.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26368,0.126341,-0.000202799,1.73568e-07,-5.94433e-11,-143626,38.0732], Tmin=(100,'K'), Tmax=(757.851,'K')), NASAPolynomial(coeffs=[13.8593,0.0407371,-2.19185e-05,4.38128e-09,-3.10406e-13,-145752,-29.5977], Tmin=(757.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1195.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFF) + group(CdCFF) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C](F)[CH]C(F)F(1426)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u1 p0 c0 {5,S} {7,S} {9,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-554.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,190,488,555,1236,1407,180,2250.74],'cm^-1')),
        HinderedRotor(inertia=(0.397589,'amu*angstrom^2'), symmetry=1, barrier=(9.14135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399336,'amu*angstrom^2'), symmetry=1, barrier=(9.18152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4542,0.0670926,-0.000128678,1.28357e-07,-4.81479e-11,-66597.3,23.7219], Tmin=(100,'K'), Tmax=(836.732,'K')), NASAPolynomial(coeffs=[4.34924,0.0286339,-1.55993e-05,3.09731e-09,-2.16322e-13,-66219.9,15.4209], Tmin=(836.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = 'F[CH]C(=C(F)F)C(F)(F)C[C](F)F(7394)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {3,S} {4,S} {9,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1139.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,190,488,555,1236,1407,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,180,1782.05],'cm^-1')),
        HinderedRotor(inertia=(0.177976,'amu*angstrom^2'), symmetry=1, barrier=(4.09203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175048,'amu*angstrom^2'), symmetry=1, barrier=(4.02469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176794,'amu*angstrom^2'), symmetry=1, barrier=(4.06485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77276,'amu*angstrom^2'), symmetry=1, barrier=(40.7592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65736,0.137615,-0.000231174,2.05369e-07,-7.18847e-11,-136834,40.9012], Tmin=(100,'K'), Tmax=(797.45,'K')), NASAPolynomial(coeffs=[13.226,0.0453433,-2.44743e-05,4.86607e-09,-3.42328e-13,-138648,-24.0219], Tmin=(797.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1139.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'FCC(=C(F)F)C(F)(F)[CH][C](F)F(7395)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {13,S}
5  F u0 p3 c0 {13,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u1 p0 c0 {8,S} {13,S} {16,S}
12 C u0 p0 c0 {6,S} {7,S} {10,D}
13 C u1 p0 c0 {4,S} {5,S} {11,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1090.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,3025,407.5,1350,352.5,182,240,577,636,1210,1413,190,488,555,1236,1407,180,180,1702.91],'cm^-1')),
        HinderedRotor(inertia=(0.203975,'amu*angstrom^2'), symmetry=1, barrier=(4.6898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204031,'amu*angstrom^2'), symmetry=1, barrier=(4.69107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201913,'amu*angstrom^2'), symmetry=1, barrier=(4.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20363,'amu*angstrom^2'), symmetry=1, barrier=(4.68185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61006,0.140956,-0.000253816,2.38184e-07,-8.61434e-11,-130974,42.1327], Tmin=(100,'K'), Tmax=(822.272,'K')), NASAPolynomial(coeffs=[10.1256,0.0504623,-2.77978e-05,5.53127e-09,-3.87532e-13,-131775,-5.32175], Tmin=(822.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1090.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C([C](F)C(F)C(F)F)=C(F)F(7396)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {15,S}
10 C u0 p0 c0 {11,S} {12,S} {13,D}
11 C u1 p0 c0 {4,S} {8,S} {10,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1152.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22178,0.124297,-0.000170934,1.22654e-07,-3.54026e-11,-138393,38.8358], Tmin=(100,'K'), Tmax=(843.229,'K')), NASAPolynomial(coeffs=[16.5238,0.0401211,-2.12002e-05,4.27695e-09,-3.07625e-13,-141385,-43.751], Tmin=(843.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1152.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CsCdCsF1s) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[C](F)C(=C(F)[CH]C(F)F)C(F)F(7397)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {15,S}
9  C u0 p0 c0 {1,S} {2,S} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {12,D} {13,S}
11 C u1 p0 c0 {9,S} {12,S} {16,S}
12 C u0 p0 c0 {5,S} {10,D} {11,S}
13 C u1 p0 c0 {6,S} {7,S} {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1185.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.976546,0.119981,-0.000170398,1.31455e-07,-4.13971e-11,-142422,38.784], Tmin=(100,'K'), Tmax=(771.14,'K')), NASAPolynomial(coeffs=[13.5337,0.0447076,-2.39657e-05,4.85047e-09,-3.48829e-13,-144660,-27.4473], Tmin=(771.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1185.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CsCFFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C(F)(F)C(F)[CH]F(7398)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
9  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {12,S} {13,D}
11 C u1 p0 c0 {4,S} {9,S} {15,S}
12 C u1 p0 c0 {5,S} {10,S} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {10,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1103.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,259,529,569,1128,1321,1390,3140,350,440,435,1725,334,575,1197,1424,3202,234,589,736,816,1240,3237,182,240,577,636,1210,1413,251.454,251.53,2038.92],'cm^-1')),
        HinderedRotor(inertia=(0.15114,'amu*angstrom^2'), symmetry=1, barrier=(6.78573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151143,'amu*angstrom^2'), symmetry=1, barrier=(6.78529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.795328,'amu*angstrom^2'), symmetry=1, barrier=(35.7149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.795671,'amu*angstrom^2'), symmetry=1, barrier=(35.7135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84157,0.139825,-0.000227118,1.93232e-07,-6.53293e-11,-132556,40.6794], Tmin=(100,'K'), Tmax=(769.077,'K')), NASAPolynomial(coeffs=[15.8048,0.0414149,-2.2248e-05,4.43138e-09,-3.12853e-13,-135074,-38.5459], Tmin=(769.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1103.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'FC=CC(F)(F)[C]([C](F)F)C(F)F(7201)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u1 p0 c0 {8,S} {9,S} {12,S}
11 C u0 p0 c0 {8,S} {13,D} {15,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 C u0 p0 c0 {7,S} {11,D} {16,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1128.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,522,611,926,1093,1137,1374,1416,3112,360,370,350,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,194,682,905,1196,1383,3221,180,180,2285.45],'cm^-1')),
        HinderedRotor(inertia=(0.403591,'amu*angstrom^2'), symmetry=1, barrier=(9.27936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402913,'amu*angstrom^2'), symmetry=1, barrier=(9.26376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403833,'amu*angstrom^2'), symmetry=1, barrier=(9.28492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403668,'amu*angstrom^2'), symmetry=1, barrier=(9.28112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75168,0.145142,-0.000265486,2.48366e-07,-8.86878e-11,-135570,41.5973], Tmin=(100,'K'), Tmax=(839.289,'K')), NASAPolynomial(coeffs=[10.5122,0.0490877,-2.66049e-05,5.23082e-09,-3.62814e-13,-136304,-7.52966], Tmin=(839.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1128.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFF) + group(CsCsFFH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]C(=CF)C(F)(F)C(F)C(F)F(7399)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
11 C u0 p0 c0 {9,S} {12,D} {13,S}
12 C u0 p0 c0 {6,S} {11,D} {16,S}
13 C u2 p0 c0 {7,S} {11,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1086.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,274,345,380,539,705,1166,1213,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,194,682,905,1196,1383,3221,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8433,0.140565,-0.000232789,2.01366e-07,-6.89126e-11,-130477,38.4995], Tmin=(100,'K'), Tmax=(782.143,'K')), NASAPolynomial(coeffs=[15.3305,0.0418724,-2.26829e-05,4.52341e-09,-3.19102e-13,-132831,-38.0104], Tmin=(782.143,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1086.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C]=C(C(F)F)C(F)(F)[CH]C(F)F(7400)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
9  C u0 p0 c0 {5,S} {6,S} {11,S} {15,S}
10 C u0 p0 c0 {3,S} {4,S} {12,S} {14,S}
11 C u0 p0 c0 {8,S} {9,S} {13,D}
12 C u1 p0 c0 {8,S} {10,S} {16,S}
13 C u1 p0 c0 {7,S} {11,D}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1058.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,171,205,598,1104,1143,1317,1411,3153,522,611,926,1093,1137,1374,1416,3112,350,440,435,1725,3025,407.5,1350,352.5,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.218435,'amu*angstrom^2'), symmetry=1, barrier=(5.02224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218158,'amu*angstrom^2'), symmetry=1, barrier=(5.01589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218169,'amu*angstrom^2'), symmetry=1, barrier=(5.01614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218367,'amu*angstrom^2'), symmetry=1, barrier=(5.02069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75883,0.142634,-0.000251973,2.31258e-07,-8.23782e-11,-127062,41.9028], Tmin=(100,'K'), Tmax=(814.779,'K')), NASAPolynomial(coeffs=[12.0829,0.0473076,-2.60857e-05,5.20052e-09,-3.65123e-13,-128409,-16.4644], Tmin=(814.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1058.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Cs_S) + radical(Cdj(Cd-CsCs)(F1s))"""),
)

species(
    label = 'FC=[C]C(F)(F)C(F)(F)[CH]C(F)F(6745)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {6,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {10,S} {15,S}
12 C u0 p0 c0 {7,S} {13,D} {16,S}
13 C u1 p0 c0 {9,S} {12,D}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1039.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,136,307,446,511,682,757,1180,1185,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,615,860,1140,1343,3152,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.161575,'amu*angstrom^2'), symmetry=1, barrier=(3.71493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16287,'amu*angstrom^2'), symmetry=1, barrier=(3.74471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162341,'amu*angstrom^2'), symmetry=1, barrier=(3.73255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452009,'amu*angstrom^2'), symmetry=1, barrier=(10.3926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3331.56,'J/mol'), sigma=(5.85389,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=520.38 K, Pc=37.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69556,0.138383,-0.000233332,2.06687e-07,-7.20314e-11,-124810,42.5185], Tmin=(100,'K'), Tmax=(798.461,'K')), NASAPolynomial(coeffs=[13.6759,0.0441309,-2.38702e-05,4.74727e-09,-3.3389e-13,-126715,-24.7378], Tmin=(798.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1039.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'F[C]F(138)',
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
    label = 'FC=[C]C(F)(F)[CH]C(F)F(3092)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-616.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.202188,'amu*angstrom^2'), symmetry=1, barrier=(4.6487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20259,'amu*angstrom^2'), symmetry=1, barrier=(4.65795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202469,'amu*angstrom^2'), symmetry=1, barrier=(4.65516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0503134,0.0992249,-0.000165008,1.49658e-07,-5.38895e-11,-73980.5,33.3083], Tmin=(100,'K'), Tmax=(785.371,'K')), NASAPolynomial(coeffs=[9.29462,0.0375631,-2.03717e-05,4.07718e-09,-2.88517e-13,-75014.6,-6.756], Tmin=(785.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-616.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCFH) + radical(Cs_S) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'FC=C1C(F)(F)C(C(F)F)C1(F)F(7193)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {8,S} {15,S}
12 C u0 p0 c0 {9,S} {10,S} {13,D}
13 C u0 p0 c0 {7,S} {12,D} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1374.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.656603,0.113129,-0.000147846,1.0496e-07,-3.07658e-11,-165108,31.1013], Tmin=(100,'K'), Tmax=(821.55,'K')), NASAPolynomial(coeffs=[13.1925,0.0457014,-2.47382e-05,5.06423e-09,-3.6782e-13,-167384,-32.9907], Tmin=(821.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1374.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCFF) + group(CsCCFF) + group(CsCsFFH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane)"""),
)

species(
    label = 'FC=C(C(F)F)C(F)(F)C=C(F)F(6752)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u0 p0 c0 {8,S} {13,D} {15,S}
12 C u0 p0 c0 {5,S} {10,D} {16,S}
13 C u0 p0 c0 {6,S} {7,S} {11,D}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1358.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42489,0.13108,-0.000213629,1.86975e-07,-6.51587e-11,-163250,37.6045], Tmin=(100,'K'), Tmax=(779.671,'K')), NASAPolynomial(coeffs=[13.0826,0.0445351,-2.38154e-05,4.74167e-09,-3.3474e-13,-165144,-26.4132], Tmin=(779.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1358.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFF)"""),
)

species(
    label = 'F[C](F)[C]1C(F)C(C(F)F)C1(F)F(7401)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {11,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {8,S} {12,S} {15,S}
11 C u0 p0 c0 {4,S} {5,S} {8,S} {16,S}
12 C u1 p0 c0 {9,S} {10,S} {13,S}
13 C u1 p0 c0 {6,S} {7,S} {12,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
"""),
    E0 = (-1140.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.79782,0.11554,-0.000167888,1.37275e-07,-4.61624e-11,-136994,36.9984], Tmin=(100,'K'), Tmax=(722.426,'K')), NASAPolynomial(coeffs=[11.6753,0.0464705,-2.44617e-05,4.90543e-09,-3.50417e-13,-138796,-19.1204], Tmin=(722.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1140.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsFFH) + group(CsCsFFH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s))"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[C]C(=C(F)F)C(F)(F)CC(F)F(7402)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {13,S}
8  C u0 p0 c0 {9,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {16,S}
11 C u0 p0 c0 {9,S} {12,D} {13,S}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 C u2 p0 c0 {7,S} {11,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
"""),
    E0 = (-1115.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,274,345,380,539,705,1166,1213,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78084,0.140986,-0.00024073,2.15083e-07,-7.5425e-11,-133912,38.4467], Tmin=(100,'K'), Tmax=(799.729,'K')), NASAPolynomial(coeffs=[13.5876,0.0449336,-2.45868e-05,4.90662e-09,-3.45549e-13,-135757,-28.427], Tmin=(799.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1115.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(CsCCl_triplet)"""),
)

species(
    label = 'FC=C(C(F)F)C(F)(F)[CH][C](F)F(7403)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {13,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u1 p0 c0 {8,S} {13,S} {15,S}
12 C u0 p0 c0 {7,S} {10,D} {16,S}
13 C u1 p0 c0 {5,S} {6,S} {11,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1107.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,190,488,555,1236,1407,180,180,1792.31],'cm^-1')),
        HinderedRotor(inertia=(0.201963,'amu*angstrom^2'), symmetry=1, barrier=(4.64353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202,'amu*angstrom^2'), symmetry=1, barrier=(4.64439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201983,'amu*angstrom^2'), symmetry=1, barrier=(4.64399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201992,'amu*angstrom^2'), symmetry=1, barrier=(4.6442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63065,0.141173,-0.000253544,2.37182e-07,-8.55984e-11,-133070,42.3141], Tmin=(100,'K'), Tmax=(821.3,'K')), NASAPolynomial(coeffs=[10.4146,0.0499666,-2.75326e-05,5.48021e-09,-3.84073e-13,-133951,-6.74411], Tmin=(821.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1107.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][CH]C(F)(F)C(=CF)C(F)(F)F(7404)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u1 p0 c0 {8,S} {13,S} {14,S}
12 C u0 p0 c0 {7,S} {10,D} {16,S}
13 C u1 p0 c0 {6,S} {11,S} {15,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {13,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1143.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,219,296,586,564,718,793,1177,1228,350,440,435,1725,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,334,575,1197,1424,3202,180,1500.24,1500.53],'cm^-1')),
        HinderedRotor(inertia=(0.285278,'amu*angstrom^2'), symmetry=1, barrier=(6.5591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28573,'amu*angstrom^2'), symmetry=1, barrier=(6.56948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285669,'amu*angstrom^2'), symmetry=1, barrier=(6.5681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286122,'amu*angstrom^2'), symmetry=1, barrier=(6.57852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6717,0.139142,-0.000239872,2.16299e-07,-7.61429e-11,-137349,41.0611], Tmin=(100,'K'), Tmax=(811.073,'K')), NASAPolynomial(coeffs=[12.7162,0.0454808,-2.46656e-05,4.89622e-09,-3.43368e-13,-138936,-20.7358], Tmin=(811.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1143.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFHH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(CdCFH) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH]C(=C(F)F)C(F)(F)C(F)C(F)F(7405)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  F u0 p3 c0 {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {15,S}
11 C u0 p0 c0 {9,S} {12,D} {13,S}
12 C u0 p0 c0 {6,S} {7,S} {11,D}
13 C u2 p0 c0 {11,S} {16,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1111.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,274,345,380,539,705,1166,1213,235,523,627,1123,1142,1372,1406,3097,350,440,435,1725,182,240,577,636,1210,1413,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65973,0.134518,-0.000194933,1.52235e-07,-4.8007e-11,-133538,38.71], Tmin=(100,'K'), Tmax=(773.771,'K')), NASAPolynomial(coeffs=[15.4355,0.0461503,-2.36376e-05,4.65909e-09,-3.29194e-13,-136184,-39.3811], Tmin=(773.771,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1111.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFF) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C(F)(F)F)C(F)(F)[CH]C(F)F(7406)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  F u0 p3 c0 {9,S}
8  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
9  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {12,S} {14,S}
11 C u0 p0 c0 {8,S} {9,S} {13,D}
12 C u1 p0 c0 {8,S} {10,S} {15,S}
13 C u1 p0 c0 {11,D} {16,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1118.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,219,296,586,564,718,793,1177,1228,522,611,926,1093,1137,1374,1416,3112,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.272384,'amu*angstrom^2'), symmetry=1, barrier=(6.26263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272276,'amu*angstrom^2'), symmetry=1, barrier=(6.26017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272406,'amu*angstrom^2'), symmetry=1, barrier=(6.26316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272306,'amu*angstrom^2'), symmetry=1, barrier=(6.26084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (208.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68896,0.137855,-0.00023116,2.0381e-07,-7.09583e-11,-134307,40.3193], Tmin=(100,'K'), Tmax=(788.625,'K')), NASAPolynomial(coeffs=[13.9984,0.0435401,-2.37196e-05,4.73875e-09,-3.34457e-13,-136323,-28.7304], Tmin=(788.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1118.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(CsCCFF) + group(CsCsFFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cds_P)"""),
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
    E0 = (-348.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-167.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-130.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-84.1182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (44.2263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (202.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-340.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-323.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-117.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-223.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-257.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-74.3385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-193.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-72.9335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-185.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-122.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (41.4532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (126.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-162.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-251.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-170.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-141.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-145.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-142.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-76.8196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-75.5657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-116.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-340.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-323.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-223.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-23.4397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-274.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-278.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-149.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-104.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-107.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['FC=C=C(F)F(1325)', 'FC(F)=CC(F)F(344)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(=C(F)F)C([C](F)F)C(F)F(6742)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC(F)=[C]C(F)C(F)(F)[CH]C(F)F(6743)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['F[CH]C(=C(F)[CH]C(F)F)C(F)(F)F(7387)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.45932e+11,'s^-1'), n=0.63878, Ea=(264.659,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(F)F-2(957)', 'F[CH]C([C](F)F)=C(F)F(3166)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(804)', 'FC(F)=[C]C(F)(F)[CH]C(F)F(7388)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['FC(F)=C1C(F)C(C(F)F)C1(F)F(7271)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['FCC(=C(F)F)C(F)(F)C=C(F)F(6751)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['FC(F)[CH]C(F)(F)[C]1C(F)C1(F)F(7389)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['F[CH][C]1C(F)(F)C(C(F)F)C1(F)F(7390)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['F[CH]C1([C](F)F)C(C(F)F)C1(F)F(7317)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(91.5016,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 91.0 to 91.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[CH]C(C(F)=CC(F)F)=C(F)F(7391)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(52.8674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C]=C(F)F(2842)', 'FC(F)=CC(F)F(344)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(3.1979,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'F[CH]C(=C(F)F)C(F)(F)C=CF(7392)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(52.2966,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'F[CH]C(=C(F)F)C(F)(F)C=C(F)F(7393)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(1.84061,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction16',
    reactants = ['FC=C=C(F)F(1325)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.13223e-11,'m^3/(mol*s)'), n=4.48095, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6004730311185978, var=1.5705211473983438, Tref=1000.0, N=276, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[C](F)[CH]C(F)F(1426)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -14.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHF(40)', 'FC(F)=[C]C(F)(F)[CH]C(F)F(7388)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)C[C](F)F(7394)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.76836e+09,'s^-1'), n=1.1815, Ea=(180.499,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FCC(=C(F)F)C(F)(F)[CH][C](F)F(7395)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(112.053,'s^-1'), n=2.49718, Ea=(42.3838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['F[CH]C([C](F)C(F)C(F)F)=C(F)F(7396)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(178.723,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['F[C](F)C(=C(F)[CH]C(F)F)C(F)F(7397)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(207.476,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)C(F)[CH]F(7398)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(161.686,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FC=CC(F)(F)[C]([C](F)F)C(F)F(7201)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(189.489,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C]C(=CF)C(F)(F)C(F)C(F)F(7399)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(213.166,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]=C(C(F)F)C(F)(F)[CH]C(F)F(7400)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C(F)(F)[CH]C(F)F(6745)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[C]F(138)', 'FC=[C]C(F)(F)[CH]C(F)F(3092)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['FC=C1C(F)(F)C(C(F)F)C1(F)F(7193)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['FC=C(C(F)F)C(F)(F)C=C(F)F(6752)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    products = ['F[C](F)[C]1C(F)C(C(F)F)C1(F)F(7401)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CF2(43)', 'FC=[C]C(F)(F)[CH]C(F)F(3092)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[C]C(=C(F)F)C(F)(F)CC(F)F(7402)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['FC=C(C(F)F)C(F)(F)[CH][C](F)F(7403)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(10500,'s^-1'), n=2.14, Ea=(33.3465,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_noH] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F[CH][CH]C(F)(F)C(=CF)C(F)(F)F(7404)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.94559e+07,'s^-1'), n=1.15307, Ea=(197.079,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]C(=C(F)F)C(F)(F)C(F)C(F)F(7405)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(211.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C(C(F)(F)F)C(F)(F)[CH]C(F)F(7406)'],
    products = ['F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #2057',
    isomers = [
        'F[CH]C(=C(F)F)C(F)(F)[CH]C(F)F(6741)',
    ],
    reactants = [
        ('FC=C=C(F)F(1325)', 'FC(F)=CC(F)F(344)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2057',
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

