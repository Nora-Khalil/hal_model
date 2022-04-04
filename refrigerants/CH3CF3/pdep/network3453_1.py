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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3442.62,'J/mol'), sigma=(5.94484,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.73 K, Pc=37.18 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.298351,0.10345,-0.000149314,1.18255e-07,-3.82483e-11,-84799.3,35.8616], Tmin=(100,'K'), Tmax=(751.605,'K')), NASAPolynomial(coeffs=[11.8112,0.0390017,-2.06903e-05,4.16338e-09,-2.9814e-13,-86619.6,-19.102], Tmin=(751.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-706.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'F[CH]C1=C[C](F)C(F)(F)C1(F)F(11916)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {13,S}
12 C u1 p0 c0 {6,S} {9,S} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-921.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0916304,0.0922566,-9.83545e-05,5.21564e-08,-1.10899e-11,-110728,29.9253], Tmin=(100,'K'), Tmax=(1128.84,'K')), NASAPolynomial(coeffs=[17.5502,0.0297429,-1.52854e-05,3.09724e-09,-2.24821e-13,-114711,-57.324], Tmin=(1128.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-921.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCCFF) + group(CsCCFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CsCdCsF1s) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2)"""),
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
    label = 'F[C]1C([C]=C(F)F)C1(F)F(11917)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-416.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,562,600,623,1070,1265,1685,370,219.052,220.74,222.425,1662.61,1665.75,1674.27],'cm^-1')),
        HinderedRotor(inertia=(0.00347825,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900582,0.0735789,-9.47962e-05,6.47649e-08,-1.79564e-11,-50023,29.771], Tmin=(100,'K'), Tmax=(873.531,'K')), NASAPolynomial(coeffs=[11.2425,0.026222,-1.34758e-05,2.70181e-09,-1.94131e-13,-51829.7,-18.7242], Tmin=(873.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-416.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Cds_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FCC(=C(F)F)C1=C(F)C1(F)F(11918)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {10,S} {11,D}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u0 p0 c0 {4,S} {7,S} {9,D}
12 C u0 p0 c0 {5,S} {6,S} {10,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-895.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.54577,0.10336,-0.000128536,7.94014e-08,-1.93636e-11,-107527,31.0922], Tmin=(100,'K'), Tmax=(1000.42,'K')), NASAPolynomial(coeffs=[18.6542,0.0265911,-1.3429e-05,2.6943e-09,-1.94519e-13,-111369,-61.5441], Tmin=(1000.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-895.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F)) + ring(Cd-Cd-Cs(F)(F))"""),
)

species(
    label = 'F[C]1C([C]2C(F)C2(F)F)C1(F)F(11919)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {10,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
11 C u1 p0 c0 {7,S} {8,S} {9,S}
12 C u1 p0 c0 {6,S} {7,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-589.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.297446,0.0886376,-0.00010588,7.06483e-08,-1.96663e-11,-70798,32.8089], Tmin=(100,'K'), Tmax=(860.305,'K')), NASAPolynomial(coeffs=[10.9157,0.0392681,-1.98013e-05,3.94473e-09,-2.82684e-13,-72625,-16.8205], Tmin=(860.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-589.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsCsFH) + ring(Cs(F)(F)-Cs(C)-Cs) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH][C]1C2C(F)(F)C2(F)C1(F)F(11920)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {12,S}
12 C u1 p0 c0 {6,S} {11,S} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-756.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.647129,0.115413,-0.000198648,1.84869e-07,-6.75302e-11,-90786,30.3262], Tmin=(100,'K'), Tmax=(798.562,'K')), NASAPolynomial(coeffs=[9.02198,0.0448766,-2.46351e-05,4.93554e-09,-3.48718e-13,-91625.5,-9.73393], Tmin=(798.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsFHH) + polycyclic(s2_3_4_ane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C([CH]C(F)=C(F)F)=C(F)F(3762)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {10,S} {11,D}
8  C u1 p0 c0 {7,S} {9,S} {13,S}
9  C u0 p0 c0 {1,S} {8,S} {12,D}
10 C u1 p0 c0 {2,S} {7,S} {14,S}
11 C u0 p0 c0 {3,S} {4,S} {7,D}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-833.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,271,519,563,612,1379,234,589,736,816,1240,3237,141,223,164,316,539,615,578,694,1133,1287,1372,1454,623.431,624.144],'cm^-1')),
        HinderedRotor(inertia=(0.000433568,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000432484,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198602,'amu*angstrom^2'), symmetry=1, barrier=(54.8971,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.308854,0.0973343,-0.000115178,6.92398e-08,-1.65669e-11,-100043,34.9682], Tmin=(100,'K'), Tmax=(1015.82,'K')), NASAPolynomial(coeffs=[17.1258,0.028682,-1.38031e-05,2.7093e-09,-1.93282e-13,-103585,-49.4175], Tmin=(1015.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-833.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(CdCFF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(C=CCJC=C) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
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
    label = 'F[CH]C(=C(F)F)C1=C(F)C1(F)F(11921)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {8,S} {11,S} {12,D}
10 C u0 p0 c0 {3,S} {7,S} {8,D}
11 C u1 p0 c0 {4,S} {9,S} {13,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-754.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,323,467,575,827,1418,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62787,'amu*angstrom^2'), symmetry=1, barrier=(37.428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62842,'amu*angstrom^2'), symmetry=1, barrier=(37.4407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.631108,0.104079,-0.000132732,8.21933e-08,-1.98883e-11,-90629.2,31.1992], Tmin=(100,'K'), Tmax=(1012.88,'K')), NASAPolynomial(coeffs=[20.1271,0.0221017,-1.13301e-05,2.28749e-09,-1.65838e-13,-94834.3,-69.2126], Tmin=(1012.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-754.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(CdCFF) + ring(Cd-Cd-Cs(F)(F)) + radical(CsCdF1sH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F))"""),
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
    label = 'F[CH]C(=C(F)F)C1C(F)=C1F(11922)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {6,S} {10,S} {11,D}
8  C u0 p0 c0 {1,S} {6,S} {9,D}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {7,S} {13,S}
11 C u0 p0 c0 {4,S} {5,S} {7,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-444.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,260,386,409,525,515,635,761,893,1354,1482,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,1465.8,1466.26],'cm^-1')),
        HinderedRotor(inertia=(0.114025,'amu*angstrom^2'), symmetry=1, barrier=(2.62165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113377,'amu*angstrom^2'), symmetry=1, barrier=(2.60676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0875527,0.0943889,-0.000134051,9.95888e-08,-2.82137e-11,-53312,30.018], Tmin=(100,'K'), Tmax=(653.83,'K')), NASAPolynomial(coeffs=[11.6462,0.0335342,-1.70572e-05,3.35971e-09,-2.37269e-13,-55034.2,-22.4458], Tmin=(653.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-444.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(CdCsCdF) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
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
    label = 'F[CH]C([C]1C(F)=C1F)=C(F)F(11923)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  C u1 p0 c0 {7,S} {8,S} {9,S}
7  C u0 p0 c0 {6,S} {10,S} {11,D}
8  C u0 p0 c0 {1,S} {6,S} {9,D}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {7,S} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {7,D}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-312.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,206,336,431,607,515,611,528,696,1312,1446,234,589,736,816,1240,3237,182,240,577,636,1210,1413,1069.39,1073.04],'cm^-1')),
        HinderedRotor(inertia=(0.165615,'amu*angstrom^2'), symmetry=1, barrier=(3.80782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00476989,'amu*angstrom^2'), symmetry=1, barrier=(3.89313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.132234,0.0921923,-0.000142337,1.19466e-07,-4.01184e-11,-37439,28.0841], Tmin=(100,'K'), Tmax=(788.465,'K')), NASAPolynomial(coeffs=[10.9113,0.0316202,-1.59005e-05,3.08854e-09,-2.15279e-13,-38955.8,-20.196], Tmin=(788.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(CdCsCdF) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Allyl_T) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
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
    label = 'FCC([C]1[C](F)C1(F)F)=C(F)F(11924)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {10,S} {13,S} {14,S}
9  C u1 p0 c0 {7,S} {10,S} {11,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 C u0 p0 c0 {5,S} {6,S} {10,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-737.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0932084,0.0929915,-0.000116404,7.91369e-08,-2.20721e-11,-88532.5,33.553], Tmin=(100,'K'), Tmax=(864.898,'K')), NASAPolynomial(coeffs=[12.294,0.0365627,-1.85347e-05,3.69583e-09,-2.64811e-13,-90642.9,-23.5373], Tmin=(864.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-737.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(Allyl_T) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]1[C](F)C1C(=C(F)F)C(F)F(11925)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {14,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u1 p0 c0 {3,S} {7,S} {11,S}
11 C u1 p0 c0 {4,S} {7,S} {10,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-686.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,165,259,333,401,397,493,1395,1505,182,240,577,636,1210,1413,180,1352.93,1352.93,1352.96,1352.97,1353.01],'cm^-1')),
        HinderedRotor(inertia=(0.156806,'amu*angstrom^2'), symmetry=1, barrier=(3.60527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156798,'amu*angstrom^2'), symmetry=1, barrier=(3.60511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.305495,0.104831,-0.000167631,1.4972e-07,-5.35747e-11,-82417.2,39.6597], Tmin=(100,'K'), Tmax=(782.12,'K')), NASAPolynomial(coeffs=[9.4304,0.0420034,-2.21346e-05,4.39133e-09,-3.09756e-13,-83541.4,-2.36899], Tmin=(782.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-686.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]C(=CF)C1C(F)(F)C1(F)F(11926)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {12,S}
11 C u0 p0 c0 {5,S} {10,D} {14,S}
12 C u2 p0 c0 {6,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-690.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,186,258,302,356,410,480,483,561,545,633,1193,1235,1433,1517,350,440,435,1725,194,682,905,1196,1383,3221,402.053,402.681,403.006,403.686,404.471,979.163,1869.88,1870.45],'cm^-1')),
        HinderedRotor(inertia=(0.00104151,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237099,'amu*angstrom^2'), symmetry=1, barrier=(27.3653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.296634,0.101359,-0.000129494,8.51086e-08,-2.24898e-11,-82845.6,32.0248], Tmin=(100,'K'), Tmax=(917.805,'K')), NASAPolynomial(coeffs=[15.4836,0.032585,-1.70932e-05,3.46393e-09,-2.50657e-13,-85742.2,-42.7524], Tmin=(917.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-690.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2)"""),
)

species(
    label = 'F[C]=C(C(F)F)C1[C](F)C1(F)F(11927)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {14,S}
10 C u0 p0 c0 {7,S} {9,S} {12,D}
11 C u1 p0 c0 {5,S} {7,S} {8,S}
12 C u1 p0 c0 {6,S} {10,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-636.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,171,205,598,1104,1143,1317,1411,3153,350,440,435,1725,212,367,445,1450,167,640,1190,382.197,382.293,382.314,382.421,1562.47,1562.73],'cm^-1')),
        HinderedRotor(inertia=(0.00115285,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00575714,'amu*angstrom^2'), symmetry=1, barrier=(9.97327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.335612,0.10352,-0.000145096,1.08078e-07,-3.25159e-11,-76366.9,36.7399], Tmin=(100,'K'), Tmax=(809.16,'K')), NASAPolynomial(coeffs=[13.4897,0.0351755,-1.84005e-05,3.6922e-09,-2.64446e-13,-78604.2,-27.0314], Tmin=(809.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-636.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Cdj(Cd-CsCs)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=[C]C(F)(F)C1[C](F)C1(F)F(11628)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
10 C u1 p0 c0 {5,S} {7,S} {8,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-603.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,136,307,446,511,682,757,1180,1185,212,367,445,1450,615,860,1140,1343,3152,1685,370,180,180,180,180,898.193,898.26],'cm^-1')),
        HinderedRotor(inertia=(0.00641377,'amu*angstrom^2'), symmetry=1, barrier=(3.67521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00641543,'amu*angstrom^2'), symmetry=1, barrier=(3.67333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3367.08,'J/mol'), sigma=(5.87221,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.93 K, Pc=37.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.267584,0.102116,-0.000135876,9.57184e-08,-2.7375e-11,-72441.9,34.1724], Tmin=(100,'K'), Tmax=(847.426,'K')), NASAPolynomial(coeffs=[13.7065,0.0361562,-1.91242e-05,3.8711e-09,-2.79295e-13,-74810.3,-30.9313], Tmin=(847.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-603.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)C1=C[C](F)C(F)(F)C1F(11928)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u1 p0 c0 {4,S} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {14,S}
12 C u1 p0 c0 {5,S} {6,S} {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-898.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.203578,0.0876136,-8.87732e-05,4.5025e-08,-9.29049e-12,-107954,30.3726], Tmin=(100,'K'), Tmax=(1148.36,'K')), NASAPolynomial(coeffs=[15.9134,0.0328924,-1.72955e-05,3.52922e-09,-2.5674e-13,-111562,-47.5912], Tmin=(1148.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-898.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCCFH) + group(CsCCFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CsCdCsF1s) + radical(Csj(Cd-CsCd)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC=[C]C1[C](F)C1(F)F(11660)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-215.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,615,860,1140,1343,3152,1685,370,431.806,432.944,432.962,433.316,1697.96,1698.1],'cm^-1')),
        HinderedRotor(inertia=(0.0046784,'amu*angstrom^2'), symmetry=1, barrier=(9.55134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25685,0.0642257,-7.31779e-05,4.4422e-08,-1.10395e-11,-25855.9,27.481], Tmin=(100,'K'), Tmax=(965.706,'K')), NASAPolynomial(coeffs=[10.7034,0.0250963,-1.23974e-05,2.46136e-09,-1.7647e-13,-27680.3,-17.7633], Tmin=(965.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Cds_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C1C2C(F)(F)C2(F)C1(F)F(11783)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {11,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {6,S} {11,D} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-996.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153898,0.0902233,-9.89159e-05,5.52268e-08,-1.25493e-11,-119776,25.6177], Tmin=(100,'K'), Tmax=(1050.64,'K')), NASAPolynomial(coeffs=[15.127,0.0332183,-1.75307e-05,3.58574e-09,-2.61381e-13,-122922,-47.3586], Tmin=(1050.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-996.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC=C(C1=C(F)C1(F)F)C(F)F(11929)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
9  C u0 p0 c0 {7,S} {10,S} {11,D}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u0 p0 c0 {5,S} {7,S} {9,D}
12 C u0 p0 c0 {6,S} {10,D} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-912.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.577882,0.103725,-0.000128837,7.92144e-08,-1.91951e-11,-109622,31.3141], Tmin=(100,'K'), Tmax=(1007.43,'K')), NASAPolynomial(coeffs=[18.9978,0.0259991,-1.31067e-05,2.62949e-09,-1.89902e-13,-113567,-63.2721], Tmin=(1007.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-912.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCFFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(CdCFH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F)) + ring(Cd-Cd-Cs(F)(F))"""),
)

species(
    label = 'F[C](F)[C]1C(F)C2(F)C1C2(F)F(11930)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {7,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,S} {14,S}
11 C u1 p0 c0 {8,S} {10,S} {12,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-738.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674111,0.0900628,-7.44483e-05,-7.0705e-08,1.21739e-10,-88675.1,27.2351], Tmin=(100,'K'), Tmax=(461.049,'K')), NASAPolynomial(coeffs=[9.84251,0.0423673,-2.28905e-05,4.57187e-09,-3.22913e-13,-89859,-13.5693], Tmin=(461.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-738.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsFFH) + polycyclic(s2_3_4_ane) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC=C([C]1[C](F)C1(F)F)C(F)F(11931)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {10,S} {13,S}
9  C u1 p0 c0 {7,S} {10,S} {11,S}
10 C u0 p0 c0 {8,S} {9,S} {12,D}
11 C u1 p0 c0 {5,S} {7,S} {9,S}
12 C u0 p0 c0 {6,S} {10,D} {14,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-754.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0790866,0.0931376,-0.000115908,7.78696e-08,-2.14189e-11,-90628.7,33.7109], Tmin=(100,'K'), Tmax=(877.094,'K')), NASAPolynomial(coeffs=[12.6007,0.0360327,-1.82479e-05,3.63932e-09,-2.60877e-13,-92825.2,-25.0564], Tmin=(877.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-754.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFH) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs-Cs(F)(F)-Cs) + radical(Allyl_T) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]C(=C(F)F)C1C(F)C1(F)F(11932)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {7,S} {8,S} {14,S}
10 C u0 p0 c0 {7,S} {11,D} {12,S}
11 C u0 p0 c0 {4,S} {5,S} {10,D}
12 C u2 p0 c0 {6,S} {10,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-689.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,350,440,435,1725,182,240,577,636,1210,1413,180,180,180,829.839,829.903,829.904,829.906,829.968],'cm^-1')),
        HinderedRotor(inertia=(0.00748247,'amu*angstrom^2'), symmetry=1, barrier=(3.6559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0074792,'amu*angstrom^2'), symmetry=1, barrier=(3.65578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0946142,0.0952954,-0.000112597,6.79309e-08,-1.65026e-11,-82736.5,32.1885], Tmin=(100,'K'), Tmax=(993.617,'K')), NASAPolynomial(coeffs=[15.7597,0.0314719,-1.62488e-05,3.28722e-09,-2.38174e-13,-85887.2,-44.1983], Tmin=(993.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-689.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCCl_triplet) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC=C(C1[C](F)[C]1F)C(F)(F)F(11933)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {12,D}
10 C u1 p0 c0 {4,S} {7,S} {11,S}
11 C u1 p0 c0 {5,S} {7,S} {10,S}
12 C u0 p0 c0 {6,S} {9,D} {14,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-734.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.269987,0.102303,-0.00015459,1.30428e-07,-4.47634e-11,-88148.2,38.1363], Tmin=(100,'K'), Tmax=(745.222,'K')), NASAPolynomial(coeffs=[10.9037,0.038912,-2.01203e-05,3.98328e-09,-2.81631e-13,-89718.7,-11.8483], Tmin=(745.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-734.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]C(=C(F)F)C1C(F)(F)C1(F)F(11934)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {10,D}
12 C u2 p0 c0 {10,S} {14,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-715.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,186,258,302,356,410,480,483,561,545,633,1193,1235,1433,1517,350,440,435,1725,182,240,577,636,1210,1413,438.058,438.059,438.067,438.073,438.073,438.079,438.082,438.102],'cm^-1')),
        HinderedRotor(inertia=(0.000878546,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000878329,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.393456,0.0991006,-0.000107241,5.97934e-08,-1.34669e-11,-85895.7,33.21], Tmin=(100,'K'), Tmax=(1067.56,'K')), NASAPolynomial(coeffs=[16.7973,0.0346916,-1.67455e-05,3.283e-09,-2.33917e-13,-89566.2,-50.8502], Tmin=(1067.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(AllylJ2_triplet) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2)"""),
)

species(
    label = '[CH]=C(C1[C](F)C1(F)F)C(F)(F)F(11935)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {12,D}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 C u1 p0 c0 {10,D} {14,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-696.434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,219,296,586,564,718,793,1177,1228,350,440,435,1725,212,367,445,1450,3120,650,792.5,1650,316.461,322.449,1435.65,1438.2,1442.26],'cm^-1')),
        HinderedRotor(inertia=(0.126849,'amu*angstrom^2'), symmetry=1, barrier=(8.89201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00158069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.211981,0.0981671,-0.00012251,7.84962e-08,-2.01606e-11,-83614.6,34.9595], Tmin=(100,'K'), Tmax=(944.956,'K')), NASAPolynomial(coeffs=[15.6819,0.0308882,-1.57129e-05,3.15047e-09,-2.26891e-13,-86618.4,-40.8193], Tmin=(944.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-696.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCdFFF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Cds_P) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-203.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-58.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (11.6076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-153.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (288.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-224.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-169.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-1.81082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-107.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-111.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-178.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-39.2676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-49.8105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (143.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (22.2302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (186.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (63.7025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (212.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-75.2607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-23.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (27.6747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (25.5013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (40.1844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (12.7669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-146.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (308.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-224.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-169.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-107.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (70.9195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-32.6137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-154.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (3.85151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-1.56983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (8.08566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FC1=CC1(F)F(974)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(29.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 29.7 kJ/mol to match endothermicity of reaction."""),
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
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[CH]C1=C[C](F)C(F)(F)C1(F)F(11916)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(79.2884,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]F(804)', 'F[C]1C([C]=C(F)F)C1(F)F(11917)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FC(F)=C1C(F)C2(F)C1C2(F)F(11844)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FCC(=C(F)F)C1=C(F)C1(F)F(11918)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[C]1C([C]2C(F)C2(F)F)C1(F)F(11919)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[CH][C]1C2C(F)(F)C2(F)C1(F)F(11920)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH]C([CH]C(F)=C(F)F)=C(F)F(3762)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[CH]C1([C](F)F)C2C(F)(F)C21F(11884)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.86885e+09,'s^-1'), n=0.611527, Ea=(54.7386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC1=CC1(F)F(974)', 'F[CH][C]=C(F)F(2842)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'F[CH]C(=C(F)F)C1=C(F)C1(F)F(11921)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(2.88655,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'F[CH]C(=C(F)F)C1C(F)=C1F(11922)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(24.3867,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FC=C=C(F)F(1325)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'F[CH]C([C]1C(F)=C1F)=C(F)F(11923)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.28222,'m^3/(mol*s)'), n=1.29695, Ea=(166.799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHF(40)', 'F[C]1C([C]=C(F)F)C1(F)F(11917)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[CH]C([C]1C(F)C1(F)F)=C(F)F(3852)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.71585e+08,'s^-1'), n=1.23602, Ea=(157.766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_noH;Cs_H_out_noH] + [R2H_S_cy3;C_rad_out_single;Cs_H_out_noH] + [R2H_S_cy3;C_rad_out_noH;Cs_H_out] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FCC([C]1[C](F)C1(F)F)=C(F)F(11924)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.18083e+09,'s^-1'), n=1.04667, Ea=(209.2,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]1[C](F)C1C(=C(F)F)C(F)F(11925)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00220343,'s^-1'), n=4.50663, Ea=(223.733,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]C(=CF)C1C(F)(F)C1(F)F(11926)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00220343,'s^-1'), n=4.50663, Ea=(225.153,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]=C(C(F)F)C1[C](F)C1(F)F(11927)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(185.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC=[C]C(F)(F)C1[C](F)C1(F)F(11628)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[C](F)C1=C[C](F)C(F)(F)C1F(11928)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(86.8774,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]F(138)', 'FC=[C]C1[C](F)C1(F)F(11660)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FC=C1C2C(F)(F)C2(F)C1(F)F(11783)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FC=C(C1=C(F)C1(F)F)C(F)F(11929)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['F[C](F)[C]1C(F)C2(F)C1C2(F)F(11930)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CF2(43)', 'FC=[C]C1[C](F)C1(F)F(11660)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction31',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FC=C([C]1[C](F)C1(F)F)C(F)F(11931)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.97418e+09,'s^-1'), n=1.23333, Ea=(200.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['F[C]C(=C(F)F)C1C(F)C1(F)F(11932)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_single;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    products = ['FC=C(C1[C](F)[C]1F)C(F)(F)F(11933)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(236.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C(=C(F)F)C1C(F)(F)C1(F)F(11934)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.026e-05,'s^-1'), n=5.16802, Ea=(223.502,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C(C1[C](F)C1(F)F)C(F)(F)F(11935)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.0279241,'s^-1'), n=4.16824, Ea=(214.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #3453',
    isomers = [
        'F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)',
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
    label = 'PDepNetwork #3453',
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

