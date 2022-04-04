species(
    label = 'FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
10 C u1 p0 c0 {7,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-572.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,2950,1000,562,600,623,1070,1265,1685,370,180,1180.3,1185.88,1187.27],'cm^-1')),
        HinderedRotor(inertia=(0.175194,'amu*angstrom^2'), symmetry=1, barrier=(4.02805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174127,'amu*angstrom^2'), symmetry=1, barrier=(4.00352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.239617,0.102848,-0.000149338,1.21667e-07,-4.09483e-11,-68670,35.6444], Tmin=(100,'K'), Tmax=(720.22,'K')), NASAPolynomial(coeffs=[10.689,0.0421536,-2.29334e-05,4.66471e-09,-3.36115e-13,-70244.3,-13.4936], Tmin=(720.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1=CC1(F)F(974)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,D} {7,S}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-329.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,323,467,575,827,1418,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2932.45,'J/mol'), sigma=(5.0155,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.04 K, Pc=52.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90533,0.00560985,9.08669e-05,-1.89954e-07,1.15706e-10,-39564.3,10.133], Tmin=(10,'K'), Tmax=(565.045,'K')), NASAPolynomial(coeffs=[4.14369,0.0253568,-1.84551e-05,6.16339e-09,-7.67829e-13,-39933.5,6.09131], Tmin=(565.045,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-329.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {12,S} {14,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-573.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,212,367,445,1450,562,600,623,1070,1265,1685,370,180,1026.36,1026.38,1026.46,1026.5,1026.6],'cm^-1')),
        HinderedRotor(inertia=(0.17847,'amu*angstrom^2'), symmetry=1, barrier=(4.10338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178438,'amu*angstrom^2'), symmetry=1, barrier=(4.10263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3376.24,'J/mol'), sigma=(5.58504,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.36 K, Pc=43.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.248578,0.091269,-0.000109832,7.22836e-08,-1.98811e-11,-68821,34.1605], Tmin=(100,'K'), Tmax=(867.757,'K')), NASAPolynomial(coeffs=[11.3411,0.040137,-2.14459e-05,4.38e-09,-3.18273e-13,-70746.2,-17.7814], Tmin=(867.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,S} {12,D}
10 C u1 p0 c0 {7,S} {8,S} {13,S}
11 C u1 p0 c0 {4,S} {9,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-706.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,350,440,435,1725,2950,1000,234,589,736,816,1240,3237,182,240,577,636,1210,1413,455.201,455.31,455.374],'cm^-1')),
        HinderedRotor(inertia=(0.000809723,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000811545,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.298351,0.10345,-0.000149314,1.18255e-07,-3.82483e-11,-84799.3,35.8616], Tmin=(100,'K'), Tmax=(751.605,'K')), NASAPolynomial(coeffs=[11.8112,0.0390017,-2.06903e-05,4.16338e-09,-2.9814e-13,-86619.6,-19.102], Tmin=(751.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-706.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[C]=C(F)F(1205)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u2 p0 c0 {3,D}
"""),
    E0 = (195.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([182,240,577,636,1210,1413],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28657,0.0148134,-1.42116e-05,6.41427e-09,-1.12277e-12,23593,10.1566], Tmin=(100,'K'), Tmax=(1382.13,'K')), NASAPolynomial(coeffs=[7.28631,0.0032378,-1.64877e-06,3.54597e-10,-2.66861e-14,22487.4,-10.4342], Tmin=(1382.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'F[CH]C1(F)[CH]C1(F)F(11484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-317.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,2950,1000,334,575,1197,1424,3202,541.842,542.001,542.237],'cm^-1')),
        HinderedRotor(inertia=(0.000573644,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15803,0.0639692,-7.44679e-05,4.36093e-08,-1.01166e-11,-38049.1,23.2582], Tmin=(100,'K'), Tmax=(1048.07,'K')), NASAPolynomial(coeffs=[13.2155,0.0179519,-8.60838e-06,1.71718e-09,-1.24049e-13,-40576.5,-35.4784], Tmin=(1048.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Csj(Cs-F1sCsCs)(F1s)(H)_1901_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC(F)=C1C(F)C2(F)C1C2(F)F(11844)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,S} {14,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-969.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375959,0.0861117,-9.02025e-05,4.85341e-08,-1.07616e-11,-116433,24.6467], Tmin=(100,'K'), Tmax=(1066.14,'K')), NASAPolynomial(coeffs=[13.9183,0.0353025,-1.8716e-05,3.83259e-09,-2.79402e-13,-119321,-41.5542], Tmin=(1066.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-969.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC(F)=C=C(F)C1(F)CC1(F)F(11845)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {7,S} {12,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {10,D} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-889.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.25784,0.101551,-0.000137342,9.8399e-08,-2.85242e-11,-106852,31.8254], Tmin=(100,'K'), Tmax=(837.941,'K')), NASAPolynomial(coeffs=[13.6743,0.0350442,-1.82878e-05,3.67913e-09,-2.64384e-13,-109187,-32.9256], Tmin=(837.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-889.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'F[C](F)C=C(F)C(F)[C]=C(F)F(3743)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {10,S} {14,S}
10 C u1 p0 c0 {3,S} {4,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-738.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,3010,987.5,1337.5,450,1655,161,297,490,584,780,1358,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.415821,'amu*angstrom^2'), symmetry=1, barrier=(9.56054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416345,'amu*angstrom^2'), symmetry=1, barrier=(9.57259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81687,'amu*angstrom^2'), symmetry=1, barrier=(41.7735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.899213,0.120467,-0.000207709,1.87771e-07,-6.6136e-11,-88621.1,36.8855], Tmin=(100,'K'), Tmax=(816.894,'K')), NASAPolynomial(coeffs=[11.1778,0.0402246,-2.16118e-05,4.26946e-09,-2.98453e-13,-89890.1,-14.6264], Tmin=(816.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-738.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cds_S)"""),
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
    label = 'FC(F)=[C]C(F)C1=CC1(F)F(11846)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-469.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,1000,562,600,623,1070,1265,1685,370,180,180,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.74428,'amu*angstrom^2'), symmetry=1, barrier=(40.1045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74659,'amu*angstrom^2'), symmetry=1, barrier=(40.1575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121189,0.0930151,-0.000131672,1.01754e-07,-3.20455e-11,-56341.4,31.6898], Tmin=(100,'K'), Tmax=(771.866,'K')), NASAPolynomial(coeffs=[11.4035,0.0345466,-1.80462e-05,3.61284e-09,-2.58172e-13,-58083.1,-19.8192], Tmin=(771.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cd-Cd(C)-Cs(F)(F)) + radical(Cds_S)"""),
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
    label = 'FC(F)=[C]C(F)C1(F)C=C1F(11847)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {6,S} {8,D}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-385.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,2950,1000,323,467,575,827,1418,562,600,623,1070,1265,1685,370,180,180,997.131],'cm^-1')),
        HinderedRotor(inertia=(0.16567,'amu*angstrom^2'), symmetry=1, barrier=(3.80908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16678,'amu*angstrom^2'), symmetry=1, barrier=(3.83459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0656608,0.0961754,-0.000154283,1.40301e-07,-5.17739e-11,-46286.7,31.8301], Tmin=(100,'K'), Tmax=(748.098,'K')), NASAPolynomial(coeffs=[8.56679,0.0405436,-2.23304e-05,4.52696e-09,-3.23898e-13,-47273.8,-4.81241], Tmin=(748.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'FC(F)=C=CC1(F)[CH]C1(F)F(11848)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {6,S} {11,D} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {9,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-472.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,180,958.244,960.365],'cm^-1')),
        HinderedRotor(inertia=(0.297099,'amu*angstrom^2'), symmetry=1, barrier=(6.8309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28672,0.0893573,-0.000124145,9.55416e-08,-3.01923e-11,-56733.9,31.8573], Tmin=(100,'K'), Tmax=(767.21,'K')), NASAPolynomial(coeffs=[10.652,0.0353119,-1.84711e-05,3.70889e-09,-2.65731e-13,-58324.3,-15.4015], Tmin=(767.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCddFF) + group(Cdd-CdsCds) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C1(F)F(11362)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-103.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2950,1000,212,367,445,1450,2726.48,2726.56],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00781,0.0264579,-1.29283e-05,-2.29098e-08,2.96552e-11,-12413.9,18.9514], Tmin=(100,'K'), Tmax=(451.813,'K')), NASAPolynomial(coeffs=[4.42599,0.0197907,-1.03421e-05,2.11919e-09,-1.54582e-13,-12602.1,12.5711], Tmin=(451.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
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
    label = 'FC(F)=C=C(F)C1(F)[CH]C1(F)F(11849)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {4,S} {7,S} {12,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {10,D} {11,D}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-663.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,145,326,398,834,1303,94,120,354,641,825,1294,540,610,2055,180,1474.26,1474.4],'cm^-1')),
        HinderedRotor(inertia=(0.198345,'amu*angstrom^2'), symmetry=1, barrier=(4.56035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.40583,0.106229,-0.000172033,1.49647e-07,-5.20051e-11,-79677.9,34.6785], Tmin=(100,'K'), Tmax=(769.565,'K')), NASAPolynomial(coeffs=[11.6119,0.0358026,-1.92421e-05,3.84301e-09,-2.71978e-13,-81291.8,-18.6204], Tmin=(769.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-663.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC#CC(F)C1(F)[CH]C1(F)F(11850)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u0 p0 c0 {4,S} {6,S} {10,S} {12,S}
9  C u1 p0 c0 {6,S} {7,S} {13,S}
10 C u0 p0 c0 {8,S} {11,T}
11 C u0 p0 c0 {5,S} {10,T}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-376.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,233,378,609,1068,1270,1314,3037,2950,1000,2175,525,239,401,1367,180,873.032,874.253,877.355],'cm^-1')),
        HinderedRotor(inertia=(0.0964506,'amu*angstrom^2'), symmetry=1, barrier=(2.21759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0996221,'amu*angstrom^2'), symmetry=1, barrier=(2.29051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.246302,0.101921,-0.000161224,1.38269e-07,-4.7385e-11,-45104,31.2686], Tmin=(100,'K'), Tmax=(782.95,'K')), NASAPolynomial(coeffs=[11.2597,0.0352118,-1.82353e-05,3.5864e-09,-2.5165e-13,-46662.8,-19.8744], Tmin=(782.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Ct-CtCs) + group(CtCF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC(F)=[C]C(F)[C]1C=C1F(11851)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
6  C u1 p0 c0 {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-45.2788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,271,519,563,612,1379,2950,1000,562,600,623,1070,1265,1685,370,180,180,1631.76,1631.9,1632.05,1632.19],'cm^-1')),
        HinderedRotor(inertia=(2.3356,'amu*angstrom^2'), symmetry=1, barrier=(53.7001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198529,'amu*angstrom^2'), symmetry=1, barrier=(4.56458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.775872,0.0782342,-0.000115877,1.02534e-07,-3.7637e-11,-5336.7,26.6857], Tmin=(100,'K'), Tmax=(736.514,'K')), NASAPolynomial(coeffs=[7.19829,0.0371487,-1.95632e-05,3.91494e-09,-2.78918e-13,-6114.43,-1.19195], Tmin=(736.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.2788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Allyl_T) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s))"""),
)

species(
    label = 'FC(F)=[C]C=C1[CH]C1(F)F(11852)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u1 p0 c0 {5,S} {6,S} {11,S}
8  C u0 p0 c0 {6,D} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-170.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2950,1000,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1685,370,512.875,512.875,512.875,512.875,512.875,512.875,512.875],'cm^-1')),
        HinderedRotor(inertia=(0.00064088,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273155,0.0712919,-6.82243e-05,3.16543e-08,-5.70253e-12,-20307.7,29.369], Tmin=(100,'K'), Tmax=(1357.8,'K')), NASAPolynomial(coeffs=[19.3996,0.0149463,-5.97747e-06,1.09156e-09,-7.52515e-14,-25501.7,-68.7551], Tmin=(1357.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(Allyl_S) + radical(Cdj(Cd-CdH)(Cd-F1sF1s))"""),
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
    label = 'FC(F)=[C]C(F)=C1[CH]C1(F)F(11853)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {7,D} {11,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-368.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2950,1000,86,203,488,582,605,741,562,600,623,1070,1265,1685,370,724.048,724.059,724.068,724.086,724.089,724.089],'cm^-1')),
        HinderedRotor(inertia=(0.000321544,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0697379,0.0767857,-7.48558e-05,3.47466e-08,-6.26148e-12,-44223.8,32.6932], Tmin=(100,'K'), Tmax=(1353.49,'K')), NASAPolynomial(coeffs=[20.7556,0.0156529,-7.1062e-06,1.37656e-09,-9.78195e-14,-49823.4,-73.3655], Tmin=(1353.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(Allyl_S) + radical(Cdj(Cd-F1sCd)(Cd-F1sF1s))"""),
)

species(
    label = 'F[C](F)C=C(F)C1(F)[CH]C1(F)F(11854)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {10,D} {12,S} {14,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-724.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.524125,0.108926,-0.000170703,1.46849e-07,-5.11406e-11,-86975.2,36.0505], Tmin=(100,'K'), Tmax=(741.061,'K')), NASAPolynomial(coeffs=[11.5121,0.0395112,-2.11955e-05,4.25137e-09,-3.02463e-13,-88637,-17.5862], Tmin=(741.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-724.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Csj(Cd-CdH)(F1s)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)[C]=C(F)C1(F)CC1(F)F(11855)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {4,S} {7,S} {12,D}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {10,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-691.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.316534,0.103937,-0.0001439,1.06326e-07,-3.19499e-11,-82985,34.3077], Tmin=(100,'K'), Tmax=(807.519,'K')), NASAPolynomial(coeffs=[13.1452,0.0372566,-2.0041e-05,4.07355e-09,-2.94097e-13,-85159.2,-27.7597], Tmin=(807.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-691.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cs) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC(F)=[C]C(F)[C]1C(F)C1(F)F(3753)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {10,S} {12,S} {14,S}
10 C u1 p0 c0 {7,S} {8,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-612.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.118758,0.102037,-0.000166143,1.55184e-07,-5.79286e-11,-73528,36.6613], Tmin=(100,'K'), Tmax=(782.672,'K')), NASAPolynomial(coeffs=[7.05673,0.0469303,-2.52013e-05,5.03813e-09,-3.56866e-13,-74086.6,7.40906], Tmin=(782.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-612.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FC(F)=C(F)C(F)[C]1[CH]C1(F)F(11856)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {13,S}
9  C u1 p0 c0 {7,S} {8,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u0 p0 c0 {4,S} {8,S} {12,D}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-620.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.715472,0.116683,-0.000200224,1.82886e-07,-6.50163e-11,-74515.8,36.5624], Tmin=(100,'K'), Tmax=(822.009,'K')), NASAPolynomial(coeffs=[9.785,0.0424193,-2.24323e-05,4.40633e-09,-3.07095e-13,-75459.4,-7.27727], Tmin=(822.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-620.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(H)_ring)"""),
)

species(
    label = 'F[C]1C(F)C1(F)C(F)[C]=C(F)F(11857)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {12,S} {14,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-547.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,259,529,569,1128,1321,1390,3140,164,312,561,654,898,1207,1299,3167,212,367,445,1450,562,600,623,1070,1265,1685,370,1127.24,1127.5],'cm^-1')),
        HinderedRotor(inertia=(0.109619,'amu*angstrom^2'), symmetry=1, barrier=(2.52035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110131,'amu*angstrom^2'), symmetry=1, barrier=(2.53213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00495699,0.0967318,-0.00012974,9.68471e-08,-2.99497e-11,-65712.9,36.2442], Tmin=(100,'K'), Tmax=(780.328,'K')), NASAPolynomial(coeffs=[10.9896,0.0404193,-2.14838e-05,4.35162e-09,-3.13763e-13,-67427.1,-14.025], Tmin=(780.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-547.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sH)(F1s)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C1(F)C(F)C(F)=C(F)F(11858)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
9  C u1 p0 c0 {4,S} {7,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 C u0 p0 c0 {3,S} {8,S} {12,D}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-570.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,284,328,853,1146,1135,1297,3239,212,367,445,1450,2950,1000,323,467,575,827,1418,182,240,577,636,1210,1413,180,180,284.616,1865.41],'cm^-1')),
        HinderedRotor(inertia=(0.0187345,'amu*angstrom^2'), symmetry=1, barrier=(46.1634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106971,'amu*angstrom^2'), symmetry=1, barrier=(2.45947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.908489,0.118065,-0.000191418,1.65362e-07,-5.69143e-11,-68500.9,35.7245], Tmin=(100,'K'), Tmax=(774.655,'K')), NASAPolynomial(coeffs=[12.8921,0.0382954,-2.04819e-05,4.07599e-09,-2.87704e-13,-70383.8,-25.6838], Tmin=(774.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-570.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-F1sCsCs)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sH)(H)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)[C]=CC1(F)C(F)C1(F)F(11859)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {4,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {7,S} {12,D} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {10,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-644.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.209319,0.100347,-0.000131878,9.14689e-08,-2.57158e-11,-77341.2,35.315], Tmin=(100,'K'), Tmax=(862.08,'K')), NASAPolynomial(coeffs=[13.8475,0.0351258,-1.83973e-05,3.71341e-09,-2.67651e-13,-79764.9,-30.4151], Tmin=(862.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-644.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cs) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsH)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)C(F)=CC1(F)[CH]C1(F)F(11860)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {13,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u0 p0 c0 {4,S} {10,D} {12,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-692.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.582825,0.110378,-0.000174549,1.51008e-07,-5.26356e-11,-83157.3,37.1885], Tmin=(100,'K'), Tmax=(756.981,'K')), NASAPolynomial(coeffs=[11.5335,0.0396585,-2.11463e-05,4.22259e-09,-2.99339e-13,-84799.9,-16.6253], Tmin=(756.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-692.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Csj(Cd-F1sCd)(F1s)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=[C]C(F)C1(F)C(F)C1(F)F(11861)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {4,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {5,S} {7,S} {11,S} {14,S}
11 C u1 p0 c0 {10,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-515.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([321,327,383,445,494,1500,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,164,312,561,654,898,1207,1299,3167,1685,370,167,640,1190,2012.12],'cm^-1')),
        HinderedRotor(inertia=(0.229854,'amu*angstrom^2'), symmetry=1, barrier=(5.28479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0167196,'amu*angstrom^2'), symmetry=1, barrier=(47.5705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.595138,0.110861,-0.000177087,1.54234e-07,-5.37828e-11,-61801.8,36.3945], Tmin=(100,'K'), Tmax=(775.643,'K')), NASAPolynomial(coeffs=[11.3557,0.039677,-2.09524e-05,4.1575e-09,-2.93273e-13,-63368.4,-16.3728], Tmin=(775.643,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=C(F)C(F)C1(F)[CH]C1(F)F(11862)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {11,S} {13,S}
10 C u1 p0 c0 {7,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {9,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-546.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,284,328,853,1146,1135,1297,3239,2950,1000,246,474,533,1155,167,640,1190,180,180,1267.73,1268.1,1268.42],'cm^-1')),
        HinderedRotor(inertia=(0.195126,'amu*angstrom^2'), symmetry=1, barrier=(4.48634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194864,'amu*angstrom^2'), symmetry=1, barrier=(4.48031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.891027,0.11744,-0.000189501,1.62816e-07,-5.57447e-11,-65602,35.7026], Tmin=(100,'K'), Tmax=(774.469,'K')), NASAPolynomial(coeffs=[13.0274,0.0377604,-2.00841e-05,3.98892e-09,-2.81287e-13,-67524.2,-26.3798], Tmin=(774.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-546.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Cdj(Cd-CsF1s)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-220.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-62.8425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-125.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (230.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-211.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-156.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-160.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (21.4637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-134.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (82.0392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (19.4357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-109.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-79.4486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (89.2835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (47.8512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (305.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (193.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (5.88064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-37.5795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-62.4269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-57.8069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-15.5636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-18.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (21.6287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-26.9397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-86.2281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (7.4756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-25.7856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['FC1=CC1(F)F(974)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=C(F)F(1205)', 'F[CH]C1(F)[CH]C1(F)F(11484)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['FC(F)=C1C(F)C2(F)C1C2(F)F(11844)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['FC(F)=C=C(F)C1(F)CC1(F)F(11845)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C](F)C=C(F)C(F)[C]=C(F)F(3743)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'FC(F)=[C]C(F)C1=CC1(F)F(11846)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(66.1346,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FC1=CC1(F)F(974)', 'F[CH][C]=C(F)F(2842)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(43.6654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'FC(F)=[C]C(F)C1(F)C=C1F(11847)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(43.1084,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC(F)=C=CC1(F)[CH]C1(F)F(11848)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(67.3208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC=C=C(F)F(1325)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(6.65884,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'FC(F)=C=C(F)C1(F)[CH]C1(F)F(11849)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(60.9,'m^3/(mol*s)'), n=1.64, Ea=(20.4815,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-1COS-R_6R!H-inRing',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-1COS-R_6R!H-inRing"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'FC#CC(F)C1(F)[CH]C1(F)F(11850)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(78)', 'FC(F)=[C]C(F)[C]1C=C1F(11851)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(7.16775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(78)', 'FC(F)=[C]C=C1[CH]C1(F)F(11852)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(20.6298,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'FC(F)=[C]C(F)=C1[CH]C1(F)F(11853)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(303.946,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[C](F)C=C(F)C1(F)[CH]C1(F)F(11854)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[C](F)[C]=C(F)C1(F)CC1(F)F(11855)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.2544e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['FC(F)=[C]C(F)[C]1C(F)C1(F)F(3753)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(162.354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['FC(F)=C(F)C(F)[C]1[CH]C1(F)F(11856)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(204.597,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]1C(F)C1(F)C(F)[C]=C(F)F(11857)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(176.514,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]1[CH]C1(F)C(F)C(F)=C(F)F(11858)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(240.578,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[C](F)[C]=CC1(F)C(F)C1(F)F(11859)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(193.221,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[C](F)C(F)=CC1(F)[CH]C1(F)F(11860)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(133.933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[C]=[C]C(F)C1(F)C(F)C1(F)F(11861)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(170.633,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[C]=C(F)C(F)C1(F)[CH]C1(F)F(11862)'],
    products = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3454',
    isomers = [
        'FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)',
    ],
    reactants = [
        ('FC1=CC1(F)F(974)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3454',
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

