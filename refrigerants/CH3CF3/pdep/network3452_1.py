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
    label = 'F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,S} {12,D}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 C u1 p0 c0 {4,S} {9,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-723.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,350,440,435,1725,212,367,445,1450,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,1013.68,1013.74,1013.77,1013.77],'cm^-1')),
        HinderedRotor(inertia=(0.00310939,'amu*angstrom^2'), symmetry=1, barrier=(2.26737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985899,'amu*angstrom^2'), symmetry=1, barrier=(2.26678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00352308,0.0934212,-0.000109755,6.69291e-08,-1.65379e-11,-86869.3,36.0761], Tmin=(100,'K'), Tmax=(975.1,'K')), NASAPolynomial(coeffs=[14.7341,0.0329936,-1.67986e-05,3.37488e-09,-2.43504e-13,-89742.1,-34.6191], Tmin=(975.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-723.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3444.78,'J/mol'), sigma=(5.95333,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=538.07 K, Pc=37.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.239617,0.102848,-0.000149338,1.21667e-07,-4.09483e-11,-68670,35.6444], Tmin=(100,'K'), Tmax=(720.22,'K')), NASAPolynomial(coeffs=[10.689,0.0421536,-2.29334e-05,4.66471e-09,-3.36115e-13,-70244.3,-13.4936], Tmin=(720.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(=C1[CH]C1(F)F)C(F)(F)F(11879)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {11,S}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 C u1 p0 c0 {7,S} {9,S} {13,S}
12 C u1 p0 c0 {6,S} {10,S} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-806.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.117713,0.0898707,-9.36182e-05,4.84005e-08,-9.96275e-12,-96899.9,32.629], Tmin=(100,'K'), Tmax=(1170.83,'K')), NASAPolynomial(coeffs=[18.2485,0.0271239,-1.32293e-05,2.62666e-09,-1.88785e-13,-101201,-58.8735], Tmin=(1170.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-806.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(CsCdFFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'F[CH]C1=C(F)[CH]C(F)(F)C1(F)F(11880)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u1 p0 c0 {8,S} {11,S} {13,S}
11 C u0 p0 c0 {5,S} {9,D} {10,S}
12 C u1 p0 c0 {6,S} {9,S} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-964.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0232948,0.0897447,-9.23969e-05,4.7104e-08,-9.62728e-12,-115855,27.7902], Tmin=(100,'K'), Tmax=(1171.45,'K')), NASAPolynomial(coeffs=[17.6948,0.0292447,-1.49287e-05,3.017e-09,-2.18618e-13,-120006,-60.4931], Tmin=(1171.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-964.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cyclopentene) + radical(cyclopentene-allyl) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2)"""),
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
    label = 'FC(F)=[C]C1(F)[CH]C1(F)F(11881)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-399.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,562,600,623,1070,1265,1685,370,180,180,180,1663.27],'cm^-1')),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.523761,0.0847709,-0.000139802,1.25215e-07,-4.45256e-11,-47949.9,29.8077], Tmin=(100,'K'), Tmax=(788.518,'K')), NASAPolynomial(coeffs=[8.91654,0.031125,-1.66907e-05,3.32336e-09,-2.34463e-13,-48929.3,-6.50611], Tmin=(788.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Cds_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC(F)=C1C(F)C2C(F)(F)C12F(11826)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,S} {14,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-996.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25208,0.0919242,-0.000110391,7.19095e-08,-1.96063e-11,-119709,23.099], Tmin=(100,'K'), Tmax=(873.484,'K')), NASAPolynomial(coeffs=[11.4708,0.0405496,-2.21675e-05,4.57488e-09,-3.34403e-13,-121669,-29.5073], Tmin=(873.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-996.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFF) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC1[C](C2(F)[CH]C2(F)F)C1(F)F(11882)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {4,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {5,S} {6,S} {8,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
11 C u1 p0 c0 {7,S} {8,S} {9,S}
12 C u1 p0 c0 {7,S} {10,S} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-595.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254256,0.102223,-0.000153386,1.30935e-07,-4.57e-11,-71527.7,34.3282], Tmin=(100,'K'), Tmax=(746.948,'K')), NASAPolynomial(coeffs=[10.1409,0.0415532,-2.15035e-05,4.26035e-09,-3.0144e-13,-72941.1,-11.8552], Tmin=(746.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFF) + group(Cs-CsCsHH) + ring(Cs(F)(F)-Cs(C)-Cs) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH][C]1C(F)(F)C2C(F)(F)C12F(11883)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {12,S}
12 C u1 p0 c0 {6,S} {11,S} {14,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-779.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.700454,0.117587,-0.00020618,1.94757e-07,-7.19499e-11,-93624,29.5013], Tmin=(100,'K'), Tmax=(798.146,'K')), NASAPolynomial(coeffs=[8.52087,0.0465325,-2.59571e-05,5.22826e-09,-3.70316e-13,-94304.8,-7.95047], Tmin=(798.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-779.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsFHH) + polycyclic(s2_3_4_ane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[CH]C(C(F)=C[C](F)F)=C(F)F(3759)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {10,S} {11,D}
8  C u0 p0 c0 {1,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 C u1 p0 c0 {2,S} {7,S} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {7,D}
12 C u1 p0 c0 {5,S} {6,S} {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-837.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,280,518,736,852,873,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,182,240,577,636,1210,1413,161,297,490,584,780,1358,1888.13],'cm^-1')),
        HinderedRotor(inertia=(0.517609,'amu*angstrom^2'), symmetry=1, barrier=(11.9008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66117,'amu*angstrom^2'), symmetry=1, barrier=(38.1936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66197,'amu*angstrom^2'), symmetry=1, barrier=(38.2119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.640855,0.110599,-0.000153735,1.10395e-07,-3.17975e-11,-100614,33.8781], Tmin=(100,'K'), Tmax=(845.7,'K')), NASAPolynomial(coeffs=[15.486,0.0343216,-1.84436e-05,3.74458e-09,-2.70231e-13,-103342,-41.2219], Tmin=(845.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-837.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + group(CdCFF) + radical(CsCdF1sH) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C1([C](F)F)C2C(F)(F)C21F(11884)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
9  C u0 p0 c0 {7,S} {8,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
11 C u1 p0 c0 {4,S} {7,S} {14,S}
12 C u1 p0 c0 {5,S} {6,S} {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-710.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.864828,0.117668,-0.000173712,1.31898e-07,-4.02118e-11,-85233.5,30.3542], Tmin=(100,'K'), Tmax=(799.625,'K')), NASAPolynomial(coeffs=[15.4448,0.0360779,-2.06529e-05,4.2835e-09,-3.12048e-13,-87841.7,-44.6827], Tmin=(799.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-710.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(CsCCCF) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsFHH) + group(CsCsFFH) + polycyclic(s2_3_3_ane) + radical(CsCsF1sH) + radical(CsCsF1sF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
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
    label = 'F[CH]C(C1=CC1(F)F)=C(F)F(11885)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {6,S} {7,D} {12,S}
10 C u1 p0 c0 {3,S} {8,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {8,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-577.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,1000,234,589,736,816,1240,3237,182,240,577,636,1210,1413,235.487,235.655,235.856,235.896,236.264,236.333],'cm^-1')),
        HinderedRotor(inertia=(0.98684,'amu*angstrom^2'), symmetry=1, barrier=(38.9827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989855,'amu*angstrom^2'), symmetry=1, barrier=(38.9865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.347842,0.0931867,-0.000107356,5.97106e-08,-1.2939e-11,-69264.5,29.4482], Tmin=(100,'K'), Tmax=(1129.57,'K')), NASAPolynomial(coeffs=[20.4984,0.0193665,-9.32664e-06,1.85394e-09,-1.33902e-13,-73973.9,-73.6624], Tmin=(1129.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-577.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cd-Cd-Cs(F)(F)) + radical(CsCdF1sH)"""),
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
    label = 'F[CH]C(=C(F)F)C1(F)C=C1F(11886)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {6,S} {10,S} {11,D}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {12,S}
10 C u1 p0 c0 {3,S} {7,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {7,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-502.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,323,467,575,827,1418,2950,1000,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.13467,'amu*angstrom^2'), symmetry=1, barrier=(3.09633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0264837,'amu*angstrom^2'), symmetry=1, barrier=(45.0581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.182172,0.100603,-0.000158303,1.36212e-07,-4.72758e-11,-60279,31.0397], Tmin=(100,'K'), Tmax=(752.684,'K')), NASAPolynomial(coeffs=[11.0117,0.0359436,-1.91379e-05,3.82156e-09,-2.71005e-13,-61817.6,-18.8109], Tmin=(752.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-502.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
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
    label = 'F[CH]C([C]1C=C1F)=C(F)F(11887)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  C u1 p0 c0 {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {11,S}
9  C u1 p0 c0 {2,S} {6,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {6,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-147.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,271,519,563,612,1379,2950,1000,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,1048.12,1048.25,1048.67,1048.9],'cm^-1')),
        HinderedRotor(inertia=(0.110777,'amu*angstrom^2'), symmetry=1, barrier=(2.54698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00327539,'amu*angstrom^2'), symmetry=1, barrier=(2.55495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848353,0.0746994,-9.31089e-05,6.48935e-08,-1.86729e-11,-17626.2,26.1276], Tmin=(100,'K'), Tmax=(838.224,'K')), NASAPolynomial(coeffs=[9.96329,0.0312046,-1.52782e-05,2.99463e-09,-2.12325e-13,-19154.3,-16.2386], Tmin=(838.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Allyl_T) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[CH]C([C]1C(F)C1(F)F)=C(F)F(3852)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {10,S}
10 C u0 p0 c0 {9,S} {11,S} {12,D}
11 C u1 p0 c0 {4,S} {10,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {10,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-790.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.200357,0.0886186,-9.65384e-05,5.51682e-08,-1.29207e-11,-94999,31.5693], Tmin=(100,'K'), Tmax=(1019.59,'K')), NASAPolynomial(coeffs=[13.8915,0.0349063,-1.75181e-05,3.50033e-09,-2.51901e-13,-97790.9,-34.7484], Tmin=(1019.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-790.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs(F)-Cs-Cs) + radical(Allyl_T) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C](F)C(=C1[CH]C1(F)F)C(F)F(11571)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
9  C u0 p0 c0 {7,S} {10,D} {11,S}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 C u1 p0 c0 {7,S} {9,S} {14,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-750.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,2950,1000,161,297,490,584,780,1358,180,180,1045.56,1045.57,1045.63,1045.67],'cm^-1')),
        HinderedRotor(inertia=(0.00913749,'amu*angstrom^2'), symmetry=1, barrier=(7.08882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.068139,'amu*angstrom^2'), symmetry=1, barrier=(52.8594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.170422,0.0884589,-9.31682e-05,4.99131e-08,-1.08789e-11,-90162.3,33.9734], Tmin=(100,'K'), Tmax=(1093.22,'K')), NASAPolynomial(coeffs=[15.3544,0.0329021,-1.69391e-05,3.42715e-09,-2.48359e-13,-93482.2,-40.6336], Tmin=(1093.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-750.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(Allyl_S) + radical(Csj(Cd-CsCd)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C1(F)[C](F)C1F(11888)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {7,S} {11,S} {12,D}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 C u1 p0 c0 {4,S} {9,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-693.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,350,440,435,1725,212,367,445,1450,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.0859152,'amu*angstrom^2'), symmetry=1, barrier=(1.97536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0514296,'amu*angstrom^2'), symmetry=1, barrier=(53.8955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.515832,0.109006,-0.000173033,1.50853e-07,-5.26581e-11,-83216.8,38.6319], Tmin=(100,'K'), Tmax=(781.671,'K')), NASAPolynomial(coeffs=[10.904,0.0401664,-2.0973e-05,4.14198e-09,-2.91386e-13,-84684.4,-11.6164], Tmin=(781.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-693.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C1(F)C(=C(F)F)C(F)F(11889)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u1 p0 c0 {4,S} {7,S} {11,S}
11 C u1 p0 c0 {7,S} {10,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-669.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,212,367,445,1450,2950,1000,182,240,577,636,1210,1413,180,1190.16,1190.86,1190.97],'cm^-1')),
        HinderedRotor(inertia=(0.21875,'amu*angstrom^2'), symmetry=1, barrier=(5.02949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218285,'amu*angstrom^2'), symmetry=1, barrier=(5.0188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.455039,0.11306,-0.000200747,1.9207e-07,-7.09449e-11,-80353.8,38.8982], Tmin=(100,'K'), Tmax=(820.717,'K')), NASAPolynomial(coeffs=[6.8343,0.0473925,-2.5641e-05,5.08364e-09,-3.56079e-13,-80535.1,11.3559], Tmin=(820.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-669.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(cyclopropane) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]C(=CF)C1(F)C(F)C1(F)F(11890)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {4,S} {7,S} {8,S} {13,S}
10 C u0 p0 c0 {7,S} {11,D} {12,S}
11 C u0 p0 c0 {5,S} {10,D} {14,S}
12 C u2 p0 c0 {6,S} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-656.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,350,440,435,1725,194,682,905,1196,1383,3221,327.824,327.903,329.053,1499.28],'cm^-1')),
        HinderedRotor(inertia=(0.13817,'amu*angstrom^2'), symmetry=1, barrier=(10.5979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459446,'amu*angstrom^2'), symmetry=1, barrier=(35.0982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.544775,0.110076,-0.000163517,1.25429e-07,-3.69392e-11,-78787.7,33.3673], Tmin=(100,'K'), Tmax=(651.068,'K')), NASAPolynomial(coeffs=[13.2404,0.0364381,-1.93303e-05,3.86616e-09,-2.75023e-13,-80817,-29.0218], Tmin=(651.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=C(C(F)F)C1(F)[CH]C1(F)F(11891)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,S} {13,S}
10 C u0 p0 c0 {7,S} {9,S} {12,D}
11 C u1 p0 c0 {7,S} {8,S} {14,S}
12 C u1 p0 c0 {6,S} {10,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-619.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,2950,1000,167,640,1190,180,180,915.631,915.76],'cm^-1')),
        HinderedRotor(inertia=(0.20879,'amu*angstrom^2'), symmetry=1, barrier=(4.80049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208875,'amu*angstrom^2'), symmetry=1, barrier=(4.80245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.719599,0.114861,-0.000190952,1.70139e-07,-6.00187e-11,-74293.5,36.7984], Tmin=(100,'K'), Tmax=(789.993,'K')), NASAPolynomial(coeffs=[11.276,0.0398656,-2.14829e-05,4.28065e-09,-3.01918e-13,-75743.9,-15.43], Tmin=(789.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Cdj(Cd-CsCs)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=[C]C(F)(F)C1(F)[CH]C1(F)F(11627)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {12,S}
10 C u1 p0 c0 {7,S} {8,S} {13,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-589.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,136,307,446,511,682,757,1180,1185,2950,1000,615,860,1140,1343,3152,1685,370,180,180,1015.8,1016.22],'cm^-1')),
        HinderedRotor(inertia=(0.164157,'amu*angstrom^2'), symmetry=1, barrier=(3.77429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163866,'amu*angstrom^2'), symmetry=1, barrier=(3.76759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3444.93,'J/mol'), sigma=(6.23988,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=538.09 K, Pc=32.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.627091,0.110966,-0.000164812,1.305e-07,-4.16892e-11,-70785.4,35.483], Tmin=(100,'K'), Tmax=(763.892,'K')), NASAPolynomial(coeffs=[13.4661,0.0371626,-1.9878e-05,4.00162e-09,-2.8623e-13,-72938.3,-28.7114], Tmin=(763.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-589.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)C1=C(F)[CH]C(F)(F)C1F(11892)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u1 p0 c0 {8,S} {11,S} {14,S}
11 C u0 p0 c0 {4,S} {9,D} {10,S}
12 C u1 p0 c0 {5,S} {6,S} {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-941.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.274151,0.0850899,-8.28266e-05,4.00334e-08,-7.86479e-12,-113081,28.2283], Tmin=(100,'K'), Tmax=(1201.63,'K')), NASAPolynomial(coeffs=[16.1516,0.032236,-1.68476e-05,3.42747e-09,-2.48755e-13,-116896,-51.2871], Tmin=(1201.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-941.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cyclopentene) + radical(cyclopentene-allyl) + radical(Csj(Cd-CsCd)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC=[C]C1(F)[CH]C1(F)F(11669)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-198.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,615,860,1140,1343,3152,1685,370,180,180,180,1530.83],'cm^-1')),
        HinderedRotor(inertia=(0.100778,'amu*angstrom^2'), symmetry=1, barrier=(2.31708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89356,0.0751717,-0.000116927,1.02726e-07,-3.64928e-11,-23783.3,27.4749], Tmin=(100,'K'), Tmax=(768.731,'K')), NASAPolynomial(coeffs=[8.12963,0.0304559,-1.5891e-05,3.15148e-09,-2.22666e-13,-24687.1,-4.17421], Tmin=(768.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Cds_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C1C(F)(F)C2C(F)(F)C12F(11718)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {6,S} {11,D} {14,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-1022.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0761036,0.0956546,-0.000118475,7.83844e-08,-2.14551e-11,-122880,24.1964], Tmin=(100,'K'), Tmax=(875.121,'K')), NASAPolynomial(coeffs=[12.4687,0.0390083,-2.13776e-05,4.41263e-09,-3.22457e-13,-125049,-33.9374], Tmin=(875.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1022.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFF) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'F[C](F)[C]1C(F)C2C(F)(F)C12F(11893)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {10,S} {12,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-763.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778372,0.0889763,-5.99957e-05,-1.18673e-07,1.69329e-10,-91692.9,25.6147], Tmin=(100,'K'), Tmax=(444.244,'K')), NASAPolynomial(coeffs=[9.78517,0.0430276,-2.35305e-05,4.70766e-09,-3.32137e-13,-92839.9,-14.4336], Tmin=(444.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsFFH) + polycyclic(s2_3_4_ane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
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
    label = 'F[C]C(=C(F)F)C1(F)CC1(F)F(11894)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
10 C u0 p0 c0 {7,S} {11,D} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 C u2 p0 c0 {6,S} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-707.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,222,329,445,522,589,1214,1475,2750,3150,900,1100,350,440,435,1725,182,240,577,636,1210,1413,347.829,347.829,347.829,347.829,347.829,1880.89,1880.89],'cm^-1')),
        HinderedRotor(inertia=(0.284787,'amu*angstrom^2'), symmetry=1, barrier=(24.45,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0155177,'amu*angstrom^2'), symmetry=1, barrier=(24.45,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.538453,0.109057,-0.000163253,1.31116e-07,-4.25774e-11,-84924.3,32.8386], Tmin=(100,'K'), Tmax=(751.264,'K')), NASAPolynomial(coeffs=[12.9288,0.0373484,-2.0068e-05,4.04721e-09,-2.89703e-13,-86947.6,-28.2808], Tmin=(751.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-707.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C(C(F)(F)F)C1(F)[CH][C]1F(11895)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u1 p0 c0 {5,S} {7,S} {11,S}
11 C u1 p0 c0 {7,S} {10,S} {13,S}
12 C u0 p0 c0 {6,S} {9,D} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-716.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.464678,0.111127,-0.000190109,1.76397e-07,-6.39207e-11,-86082.9,37.5329], Tmin=(100,'K'), Tmax=(813.833,'K')), NASAPolynomial(coeffs=[8.42054,0.0440948,-2.35014e-05,4.6449e-09,-3.25337e-13,-86755.5,1.25028], Tmin=(813.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(cyclopropane) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]C(=C(F)F)C1(F)C(F)C1(F)F(11896)',
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
10 C u0 p0 c0 {7,S} {11,D} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {10,D}
12 C u2 p0 c0 {10,S} {14,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-681.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,350,440,435,1725,182,240,577,636,1210,1413,311.294,311.294,311.294,311.295],'cm^-1')),
        HinderedRotor(inertia=(0.737349,'amu*angstrom^2'), symmetry=1, barrier=(50.7038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.737342,'amu*angstrom^2'), symmetry=1, barrier=(50.7038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.492565,0.106023,-0.000135205,9.3621e-08,-2.64627e-11,-81843.8,34.0237], Tmin=(100,'K'), Tmax=(856.351,'K')), NASAPolynomial(coeffs=[13.5907,0.0402421,-1.99865e-05,3.92675e-09,-2.78592e-13,-84255.9,-31.7363], Tmin=(856.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-681.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(AllylJ2_triplet) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]=C(C(F)(F)F)C1(F)[CH]C1(F)F(11897)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {12,D}
11 C u1 p0 c0 {7,S} {8,S} {13,S}
12 C u1 p0 c0 {10,D} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-679.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,219,296,586,564,718,793,1177,1228,350,440,435,1725,2950,1000,3120,650,792.5,1650,180,180,834.525],'cm^-1')),
        HinderedRotor(inertia=(0.140539,'amu*angstrom^2'), symmetry=1, barrier=(3.23127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145107,'amu*angstrom^2'), symmetry=1, barrier=(3.33629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.52336,0.108382,-0.000163074,1.3163e-07,-4.28685e-11,-81544.2,34.7743], Tmin=(100,'K'), Tmax=(749.678,'K')), NASAPolynomial(coeffs=[12.9391,0.0365638,-1.94012e-05,3.88886e-09,-2.77237e-13,-83563.1,-26.2979], Tmin=(749.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Cds_P) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-226.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-81.4834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-10.1975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (21.0724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-177.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (282.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-230.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-7.58561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-113.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-144.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-174.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (21.2241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-60.9463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (76.8712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-0.696206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (163.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (316.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (206.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-96.6192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-33.4331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-44.0413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (21.8999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (28.4627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (34.4097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (3.51666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-171.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (302.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-230.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-113.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (65.1447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-195.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-1.92327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1.16453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (2.31087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['FC1=CC1(F)F(974)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(12.5719,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 12.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
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
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[CH]C(=C1[CH]C1(F)F)C(F)(F)F(11879)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(259.874,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[CH]C1=C(F)[CH]C(F)(F)C1(F)F(11880)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(60.9712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]F(804)', 'FC(F)=[C]C1(F)[CH]C1(F)F(11881)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['FC(F)=C1C(F)C2C(F)(F)C12F(11826)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['FC1[C](C2(F)[CH]C2(F)F)C1(F)F(11882)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[CH][C]1C(F)(F)C2C(F)(F)C12F(11883)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C(C(F)=C[C](F)F)=C(F)F(3759)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[CH]C1([C](F)F)C2C(F)(F)C21F(11884)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(64.7047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[CH]C(C1=CC1(F)F)=C(F)F(11885)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(58.0723,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['FC1=CC1(F)F(974)', 'F[CH][C]=C(F)F(2842)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(1.24768,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'F[CH]C(=C(F)F)C1(F)C=C1F(11886)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(38.8738,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FC=C=C(F)F(1325)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(78)', 'F[CH]C([C]1C=C1F)=C(F)F(11887)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(5.07996,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHF(40)', 'FC(F)=[C]C1(F)[CH]C1(F)F(11881)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[CH]C([C]1C(F)C1(F)F)=C(F)F(3852)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(142.183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[C](F)C(=C1[CH]C1(F)F)C(F)F(11571)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(205.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]C(=C(F)F)C1(F)[C](F)C1F(11888)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(181.659,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]1[CH]C1(F)C(=C(F)F)C(F)F(11889)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(223.733,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]C(=CF)C1(F)C(F)C1(F)F(11890)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(217.345,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]=C(C(F)F)C1(F)[CH]C1(F)F(11891)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C1(F)[CH]C1(F)F(11627)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[C](F)C1=C(F)[CH]C(F)(F)C1F(11892)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(67.4215,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[C]F(138)', 'FC=[C]C1(F)[CH]C1(F)F(11669)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['FC=C1C(F)(F)C2C(F)(F)C12F(11718)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['F[C](F)[C]1C(F)C2C(F)(F)C12F(11893)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CF2(43)', 'FC=[C]C1(F)[CH]C1(F)F(11669)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[C]C(=C(F)F)C1(F)CC1(F)F(11894)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    products = ['FC=C(C(F)(F)F)C1(F)[CH][C]1F(11895)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(236.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]C(=C(F)F)C1(F)C(F)C1(F)F(11896)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(215.467,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(C(F)(F)F)C1(F)[CH]C1(F)F(11897)'],
    products = ['F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #3452',
    isomers = [
        'F[CH]C(=C(F)F)C1(F)[CH]C1(F)F(11572)',
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
    label = 'PDepNetwork #3452',
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

