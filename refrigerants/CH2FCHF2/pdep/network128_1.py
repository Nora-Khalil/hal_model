species(
    label = '[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
7  C u1 p0 c0 {5,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {11,S}
9  C u0 p0 c0 {4,S} {7,S} {8,D}
10 C u1 p0 c0 {6,S} {13,S} {14,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-406.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2750,3150,900,1100,271,519,563,612,1379,3000,3100,440,815,1455,1000,444.631,444.668,444.668,444.686,444.69,444.722],'cm^-1')),
        HinderedRotor(inertia=(0.000852524,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000852737,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.325484,0.0822596,-8.43213e-05,4.41016e-08,-9.31942e-12,-48714.4,26.441], Tmin=(100,'K'), Tmax=(1132.25,'K')), NASAPolynomial(coeffs=[15.3133,0.029312,-1.41781e-05,2.80237e-09,-2.00772e-13,-52108.5,-47.7285], Tmin=(1132.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = 'FC1=CC(F)=C1(307)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u0 p0 c0 {2,S} {3,D} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (29.0416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,180,1170.89,2106.02,2106.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89011,0.00686793,0.000107406,-2.40335e-07,1.58082e-10,3500.24,9.34061], Tmin=(10,'K'), Tmax=(520.912,'K')), NASAPolynomial(coeffs=[4.12106,0.0279909,-1.93509e-05,6.26873e-09,-7.66374e-13,3165.53,5.39524], Tmin=(520.912,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(29.0416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCC(F)DC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)(F)C1[C](F)C=C1F(591)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u1 p0 c0 {3,S} {5,S} {9,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {12,S}
10 C u1 p0 c0 {6,S} {13,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-422.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0821443,0.0819567,-7.96996e-05,3.81906e-08,-7.25134e-12,-50631,29.5473], Tmin=(100,'K'), Tmax=(1268.74,'K')), NASAPolynomial(coeffs=[18.4698,0.0239851,-1.11611e-05,2.17653e-09,-1.54907e-13,-55296.8,-63.5393], Tmin=(1268.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-422.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCCFH) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)[C](F)C1C=C1F(650)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 C u1 p0 c0 {6,S} {13,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-219.695,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,215,315,519,588,595,1205,1248,212,367,445,1450,323,467,575,827,1418,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0248967,0.0949074,-0.000133639,1.04177e-07,-3.3142e-11,-26286.9,32.9225], Tmin=(100,'K'), Tmax=(764.389,'K')), NASAPolynomial(coeffs=[11.2516,0.0361587,-1.8354e-05,3.63018e-09,-2.57556e-13,-28003.3,-18.2238], Tmin=(764.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-219.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd-Cd-Cs(C-F)) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'F[C](F)CC1(F)[CH]C(F)=C1(592)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {5,S} {9,S} {14,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {7,S} {8,D}
10 C u1 p0 c0 {3,S} {4,S} {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-403.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,271,519,563,612,1379,190,488,555,1236,1407,503.026,503.849,503.886,504.011,504.112,504.332,504.366],'cm^-1')),
        HinderedRotor(inertia=(0.00178879,'amu*angstrom^2'), symmetry=1, barrier=(0.32154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000751926,'amu*angstrom^2'), symmetry=1, barrier=(0.135526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3492.23,'J/mol'), sigma=(5.86594,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.48 K, Pc=39.26 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.766348,0.0771745,-8.01144e-05,4.55488e-08,-1.09151e-11,-48434.9,26.201], Tmin=(100,'K'), Tmax=(983.066,'K')), NASAPolynomial(coeffs=[10.7133,0.0367002,-1.83556e-05,3.66589e-09,-2.63691e-13,-50390.6,-21.6172], Tmin=(983.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(F)C1=CC(F)(F)[CH]1(651)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  C u1 p0 c0 {5,S} {7,S} {12,S}
9  C u0 p0 c0 {5,S} {7,D} {11,S}
10 C u1 p0 c0 {6,S} {13,S} {14,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-404.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([207,289,297,369,408,524,563,645,605,763,737,855,975,1147,1148,1250,2750,3150,900,1100,3000,3100,440,815,1455,1000,180,180,657.003,1426.63,1426.76,1426.81,1426.81,1426.91],'cm^-1')),
        HinderedRotor(inertia=(0.239648,'amu*angstrom^2'), symmetry=1, barrier=(5.50997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36213,'amu*angstrom^2'), symmetry=1, barrier=(31.3181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621718,0.0765516,-7.12628e-05,3.34984e-08,-6.44277e-12,-48578.2,24.9731], Tmin=(100,'K'), Tmax=(1221.28,'K')), NASAPolynomial(coeffs=[14.4563,0.0312404,-1.56114e-05,3.12e-09,-2.24297e-13,-51957.4,-44.5365], Tmin=(1221.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFF) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(cyclobutene-allyl) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
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
    label = 'F[C](F)C1(F)[CH]C(F)=C1(652)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-358.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,190,488,555,1236,1407,528.895,528.903,528.909,528.914,528.917,528.927],'cm^-1')),
        HinderedRotor(inertia=(0.000602608,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20972,0.063219,-6.43813e-05,3.34833e-08,-7.07243e-12,-42977.5,22.3914], Tmin=(100,'K'), Tmax=(1128.21,'K')), NASAPolynomial(coeffs=[12.3983,0.0235505,-1.16408e-05,2.31868e-09,-1.66698e-13,-45502.1,-32.937], Tmin=(1128.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'FC1=CC2(F)C1CC2(F)F(594)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-614.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905494,0.0712774,-6.18608e-05,2.72074e-08,-4.9598e-12,-73778.1,20.1125], Tmin=(100,'K'), Tmax=(1266.8,'K')), NASAPolynomial(coeffs=[13.077,0.0328449,-1.63532e-05,3.25852e-09,-2.33523e-13,-76861.8,-41.4865], Tmin=(1266.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-614.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)C2[C](F)C21(653)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  C u1 p0 c0 {4,S} {6,S} {7,S}
10 C u1 p0 c0 {8,S} {13,S} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-323.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.298927,0.0889739,-0.000109049,7.18794e-08,-1.95164e-11,-38780.5,26.9255], Tmin=(100,'K'), Tmax=(884.582,'K')), NASAPolynomial(coeffs=[12.0801,0.0357005,-1.87124e-05,3.79696e-09,-2.74988e-13,-40864.8,-28.4669], Tmin=(884.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + polycyclic(s2_3_3_ane) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'FC12[CH]C(F)([CH]1)C(F)(F)C2(633)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {6,S} {14,S}
10 C u1 p0 c0 {5,S} {6,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-367.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52645,0.063806,-5.33418e-05,2.58228e-08,-5.95119e-12,-44067.7,21.7151], Tmin=(100,'K'), Tmax=(934.579,'K')), NASAPolynomial(coeffs=[5.79313,0.0455443,-2.40313e-05,4.91433e-09,-3.58075e-13,-44865.2,1.41959], Tmin=(934.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsCsFF) + polycyclic(s3_4_5_ane) + radical(bicyclo[2.1.1]hexane-C5) + radical(bicyclo[2.1.1]hexane-C5) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C2(F)C1CC2(F)F(612)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {10,S}
10 C u1 p0 c0 {6,S} {9,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-348.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32761,0.0690109,-5.59595e-05,-2.8352e-09,2.88893e-11,-41830.9,23.5128], Tmin=(100,'K'), Tmax=(496.297,'K')), NASAPolynomial(coeffs=[5.69992,0.0461566,-2.4318e-05,4.94411e-09,-3.58326e-13,-42417.4,3.94535], Tmin=(496.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(Cs-CsCsHH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-CsCsHH) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(bicyclo[2.2.0]hexane-secondary) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]=C(F)C=C(F)C([CH2])(F)F(654)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {11,S}
8  C u1 p0 c0 {5,S} {12,S} {13,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u1 p0 c0 {9,D} {14,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-304.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,323,467,575,827,1418,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,250,446,589,854,899,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.577587,'amu*angstrom^2'), symmetry=1, barrier=(13.2799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577504,'amu*angstrom^2'), symmetry=1, barrier=(13.278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577366,'amu*angstrom^2'), symmetry=1, barrier=(13.2748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.525514,0.107381,-0.000160598,1.25832e-07,-3.93383e-11,-36498.5,31.0411], Tmin=(100,'K'), Tmax=(783.439,'K')), NASAPolynomial(coeffs=[14.2069,0.032162,-1.65827e-05,3.28322e-09,-2.32402e-13,-38806.9,-36.4389], Tmin=(783.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-CdsHH) + radical(Csj(Cs-F1sF1sCd)(H)(H)) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[CH2]C(F)(F)C1=CC(F)=C1(655)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {5,D} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {6,S} {7,D}
9  C u1 p0 c0 {4,S} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-61.5438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2750,3150,900,1100,280,518,736,852,873,3000,3100,440,815,1455,1000,180,180,180,180,1429.57,1429.94,1429.98,4000],'cm^-1')),
        HinderedRotor(inertia=(0.15413,'amu*angstrom^2'), symmetry=1, barrier=(3.54376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08525,'amu*angstrom^2'), symmetry=1, barrier=(24.952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.522064,0.0808109,-9.62584e-05,5.97794e-08,-1.49857e-11,-7280.34,25.1947], Tmin=(100,'K'), Tmax=(964.316,'K')), NASAPolynomial(coeffs=[13.3107,0.0277631,-1.37415e-05,2.73207e-09,-1.96058e-13,-9746.78,-36.0383], Tmin=(964.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.5438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + ring(Cd-Cd-Cd-Cd(F)) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
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
    label = 'C=C(F)C1(F)[CH]C(F)=C1(656)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u1 p0 c0 {4,S} {8,S} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {10,S}
7  C u0 p0 c0 {2,S} {4,S} {9,D}
8  C u0 p0 c0 {3,S} {5,S} {6,D}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-275.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,323,467,575,827,1418,271,519,563,612,1379,2950,3100,1380,975,1025,1650,543.812,544.226,544.761,544.954,545.655,545.73],'cm^-1')),
        HinderedRotor(inertia=(0.000565032,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07341,0.0639411,-5.4493e-05,2.3608e-08,-4.18793e-12,-33052.2,21.8956], Tmin=(100,'K'), Tmax=(1316.63,'K')), NASAPolynomial(coeffs=[13.1341,0.0273001,-1.27489e-05,2.47119e-09,-1.74491e-13,-36228.1,-39.6076], Tmin=(1316.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl)"""),
)

species(
    label = 'F[C]1[CH]C(F)=C1(615)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (86.8055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,271,519,563,612,1379,180,180,1104.81,1106.27,1801.34],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06134,0.0383117,-3.38918e-05,1.45559e-08,-2.45109e-12,10513.6,15.9062], Tmin=(100,'K'), Tmax=(1431.45,'K')), NASAPolynomial(coeffs=[12.0907,0.0102863,-4.52462e-06,8.78888e-10,-6.24417e-14,7642.3,-36.0769], Tmin=(1431.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.8055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(Cd-CdH)_ring) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring)"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)C=C=C1(657)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
6  C u0 p0 c0 {4,S} {9,D} {10,S}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u1 p0 c0 {5,S} {12,S} {13,S}
9  C u0 p0 c0 {6,D} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (6.3639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2750,3150,900,1100,3000,3100,440,815,1455,1000,180,180,180,180,180,180,1011.5,1011.76],'cm^-1')),
        HinderedRotor(inertia=(0.00655911,'amu*angstrom^2'), symmetry=1, barrier=(4.76699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0065698,'amu*angstrom^2'), symmetry=1, barrier=(4.7702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205809,0.0935556,-0.0001544,1.41602e-07,-5.14176e-11,892.237,25.617], Tmin=(100,'K'), Tmax=(797.529,'K')), NASAPolynomial(coeffs=[7.92144,0.0385944,-2.04396e-05,4.05327e-09,-2.85293e-13,178.769,-6.6183], Tmin=(797.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.3639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C(F)=C1[CH]C(F)=C1(658)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {7,D}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,D}
7  C u0 p0 c0 {2,S} {3,D} {8,S}
8  C u1 p0 c0 {7,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-13.5177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,206,336,431,607,515,611,528,696,1312,1446,3000,3100,440,815,1455,1000,764.393,764.394,764.396,764.401,764.404,764.405,764.408,764.413,764.42],'cm^-1')),
        HinderedRotor(inertia=(0.000288498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.093,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52233,0.0330812,4.90066e-05,-9.9536e-08,4.41316e-11,-1516.79,19.9426], Tmin=(100,'K'), Tmax=(926.624,'K')), NASAPolynomial(coeffs=[19.9871,0.00841481,-1.63784e-07,-5.65356e-11,-2.59235e-15,-7301.77,-80.4826], Tmin=(926.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.5177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(C=CC=CCJ)"""),
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
    label = '[CH2]C(F)(F)C1=C=C(F)[CH]1(659)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {9,D}
6  C u1 p0 c0 {5,S} {7,S} {10,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u1 p0 c0 {4,S} {11,S} {12,S}
9  C u0 p0 c0 {5,D} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (140.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2950,1000,132,421,945,1169,1268,3000,3100,440,815,1455,1000,180,180,972.444,972.444,972.467,972.473,972.518],'cm^-1')),
        HinderedRotor(inertia=(0.186515,'amu*angstrom^2'), symmetry=1, barrier=(4.28834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186484,'amu*angstrom^2'), symmetry=1, barrier=(4.28764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12193,0.0691263,-7.80239e-05,4.88131e-08,-1.28227e-11,16940.1,24.9304], Tmin=(100,'K'), Tmax=(904.971,'K')), NASAPolynomial(coeffs=[9.5673,0.0317991,-1.61561e-05,3.23863e-09,-2.33185e-13,15411.5,-14.9707], Tmin=(904.971,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CdCddCF) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(C=CCJC=C) + radical(Csj(Cs-F1sF1sCd)(H)(H))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)[C]=C=C1(660)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {9,D} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  C u1 p0 c0 {4,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (244.206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,3000,3100,440,815,1455,1000,180,180,180,180,180,1170.72,1170.95],'cm^-1')),
        HinderedRotor(inertia=(0.122324,'amu*angstrom^2'), symmetry=1, barrier=(2.81246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122239,'amu*angstrom^2'), symmetry=1, barrier=(2.81051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.135716,0.0985962,-0.000180696,1.72647e-07,-6.29804e-11,29497.2,27.0347], Tmin=(100,'K'), Tmax=(833.669,'K')), NASAPolynomial(coeffs=[6.96942,0.0376769,-2.04706e-05,4.04295e-09,-2.81634e-13,29335.3,1.17216], Tmin=(833.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)[C]=C(F)C1(661)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u1 p0 c0 {7,S} {13,S} {14,S}
10 C u1 p0 c0 {5,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-318.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,215,315,519,588,595,1205,1248,246,474,533,1155,3000,3100,440,815,1455,1000,277.748,313.02,367.582,448.836,459.667,1616.06,1619.64],'cm^-1')),
        HinderedRotor(inertia=(0.00381952,'amu*angstrom^2'), symmetry=1, barrier=(0.1252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589427,'amu*angstrom^2'), symmetry=1, barrier=(31.9418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.114169,0.0916525,-0.00011574,7.79447e-08,-2.13036e-11,-38133.9,28.7941], Tmin=(100,'K'), Tmax=(885.962,'K')), NASAPolynomial(coeffs=[13.0965,0.0330403,-1.65074e-05,3.27634e-09,-2.34259e-13,-40434.3,-32.2666], Tmin=(885.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'CC(F)(F)C1(F)[C]=C(F)[CH]1(662)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
8  C u1 p0 c0 {5,S} {9,S} {14,S}
9  C u0 p0 c0 {4,S} {8,S} {10,D}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-365.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,222,329,445,522,589,1214,1475,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,1000,271,519,563,612,1379,180,180,723.455,724.84,728.858],'cm^-1')),
        HinderedRotor(inertia=(0.126255,'amu*angstrom^2'), symmetry=1, barrier=(2.90286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76003,'amu*angstrom^2'), symmetry=1, barrier=(40.4666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493126,0.0827147,-9.18041e-05,5.4822e-08,-1.35276e-11,-43862.5,24.6907], Tmin=(100,'K'), Tmax=(967.323,'K')), NASAPolynomial(coeffs=[12.1205,0.0346333,-1.72443e-05,3.43542e-09,-2.46761e-13,-46112,-31.0181], Tmin=(967.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C(F)(F)[C]1C=C(F)C1F(663)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {12,S}
10 C u1 p0 c0 {6,S} {13,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-408.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684369,0.0728697,-6.41354e-05,2.84071e-08,-5.12973e-12,-48988.9,25.3445], Tmin=(100,'K'), Tmax=(1298.41,'K')), NASAPolynomial(coeffs=[14.6985,0.0296966,-1.42595e-05,2.79856e-09,-1.99012e-13,-52628.1,-45.9255], Tmin=(1298.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-408.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd(F)-Cd-Cs-Cs) + radical(Allyl_T) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'FCC(F)(F)[C]1[CH]C(F)=C1(664)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
7  C u1 p0 c0 {5,S} {8,S} {9,S}
8  C u1 p0 c0 {7,S} {10,S} {14,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {4,S} {8,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-439.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294804,0.0703826,-5.5846e-05,2.14836e-08,-3.26557e-12,-52750.1,24.3703], Tmin=(100,'K'), Tmax=(1568.91,'K')), NASAPolynomial(coeffs=[19.3616,0.0217706,-9.36857e-06,1.73405e-09,-1.18509e-13,-58732.9,-76.2028], Tmin=(1568.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-439.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd(F)-Cd-Cs-Cs) + radical(Allyl_T) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH2][C](F)C1(F)C=C(F)C1F(665)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {4,S} {6,S} {7,D}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u1 p0 c0 {9,S} {13,S} {14,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-326.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,2950,1000,323,467,575,827,1418,212,367,445,1450,3000,3100,440,815,1455,1000,180,180,972.878,973.057],'cm^-1')),
        HinderedRotor(inertia=(0.103168,'amu*angstrom^2'), symmetry=1, barrier=(2.37204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102214,'amu*angstrom^2'), symmetry=1, barrier=(2.3501,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582713,0.0854299,-0.000102561,5.41621e-08,-1.54696e-12,-39104.6,30.1096], Tmin=(100,'K'), Tmax=(563.476,'K')), NASAPolynomial(coeffs=[9.17407,0.0394629,-2.01823e-05,4.00732e-09,-2.84977e-13,-40311.2,-8.52647], Tmin=(563.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC[C](F)C1(F)[CH]C(F)=C1(666)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u1 p0 c0 {5,S} {10,S} {14,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {4,S} {8,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-382.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,2750,3150,900,1100,271,519,563,612,1379,421.487,421.846,422.349,422.354,423.241,423.368,1858.86],'cm^-1')),
        HinderedRotor(inertia=(0.281039,'amu*angstrom^2'), symmetry=1, barrier=(35.5676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0200986,'amu*angstrom^2'), symmetry=1, barrier=(35.6251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.742469,0.0777033,-8.07388e-05,4.58353e-08,-1.09554e-11,-45888.7,26.6011], Tmin=(100,'K'), Tmax=(985.91,'K')), NASAPolynomial(coeffs=[10.8431,0.0367228,-1.83888e-05,3.67412e-09,-2.64343e-13,-47880.4,-21.9853], Tmin=(985.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-382.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCsCsF1s) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)C=[C]C1F(667)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
7  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {10,D} {12,S}
9  C u1 p0 c0 {6,S} {13,S} {14,S}
10 C u1 p0 c0 {7,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-283.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,2950,1000,3000,3100,440,815,1455,1000,180,180,180,1022.14,1022.99],'cm^-1')),
        HinderedRotor(inertia=(0.0404579,'amu*angstrom^2'), symmetry=1, barrier=(30.0266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193278,'amu*angstrom^2'), symmetry=1, barrier=(4.44385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30203,0.0867507,-0.000100141,6.08985e-08,-1.51063e-11,-34005.6,28.3683], Tmin=(100,'K'), Tmax=(968.577,'K')), NASAPolynomial(coeffs=[13.3881,0.0327075,-1.6445e-05,3.29064e-09,-2.36947e-13,-36540.5,-34.3467], Tmin=(968.577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCC(F)(F)C1(F)[CH][C]=C1(668)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {10,S} {14,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-326.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,222,329,445,522,589,1214,1475,528,1116,1182,1331,1402,1494,3075,3110,2750,3150,900,1100,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542154,0.0791787,-7.95942e-05,4.16539e-08,-8.93914e-12,-39105,25.1193], Tmin=(100,'K'), Tmax=(1105.67,'K')), NASAPolynomial(coeffs=[13.5695,0.0320494,-1.56565e-05,3.10244e-09,-2.22388e-13,-41985.8,-39.0388], Tmin=(1105.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
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
    E0 = (-141.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-62.4999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (131.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-57.5586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (94.5393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (214.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-206.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-1.21453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-35.1742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-123.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (24.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (217.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (111.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (46.6054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-80.5706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (278.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (169.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (191.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (201.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (360.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (14.4848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-130.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-34.4425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-4.96337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (54.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-48.9003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (40.1349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (25.8566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['CH2CF2(57)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(73.5503,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 73.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['[CH2]C(F)(F)C1[C](F)C=C1F(591)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.53073e+11,'s^-1'), n=0.611, Ea=(152.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCsCJ;CsJ-CdH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CdH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)(F)[C](F)C1C=C1F(650)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['F[C](F)CC1(F)[CH]C(F)=C1(592)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)(F)C1=CC(F)(F)[CH]1(651)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(308.188,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(18)', 'F[C](F)C1(F)[CH]C(F)=C1(652)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['FC1=CC2(F)C1CC2(F)F(594)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.6e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['[CH2]C(F)(F)C1(F)C2[C](F)C21(653)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['FC12[CH]C(F)([CH]1)C(F)(F)C2(633)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.96041e+10,'s^-1'), n=0.12, Ea=(179.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['F[C]1[CH]C2(F)C1CC2(F)F(612)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.10979e+11,'s^-1'), n=0.118526, Ea=(91.833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra_cs2H] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(F)C=C(F)C([CH2])(F)F(654)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[CH2]C(F)(F)C1=CC(F)=C1(655)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(15.2165,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C](F)F(163)', 'FC1=CC(F)=C1(307)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000404,'m^3/(mol*s)'), n=2.41, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-2R!H-R_Ext-1R!H-R_Ext-3C-R_Ext-3C-R_Sp-6R!H-2R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-2R!H-R_Ext-1R!H-R_Ext-3C-R_Ext-3C-R_Sp-6R!H-2R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'C=C(F)C1(F)[CH]C(F)=C1(656)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(58.1543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2CF2(57)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00100541,'m^3/(mol*s)'), n=2.87982, Ea=(2.99211,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', '[CH2]C(F)(F)C1(F)C=C=C1(657)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(8.44811,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C](F)F(163)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(78)', '[CH2]C(F)=C1[CH]C(F)=C1(658)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(22.3403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[CH2]C(F)(F)C1=C=C(F)[CH]1(659)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(151.685,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', '[CH2]C(F)(F)C1(F)[C]=C=C1(660)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(206.5,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(F)(F)C1(F)[C]=C(F)C1(661)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(F)(F)C1(F)[C]=C(F)[CH]1(662)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['[CH2]C(F)(F)[C]1C=C(F)C1F(663)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(180.434,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    products = ['FCC(F)(F)[C]1[CH]C(F)=C1(664)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(209.914,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C](F)C1(F)C=C(F)C1F(665)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(189.394,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC[C](F)C1(F)[CH]C(F)=C1(666)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(142.326,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(F)(F)C1(F)C=[C]C1F(667)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(132.695,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['FCC(F)(F)C1(F)[CH][C]=C1(668)'],
    products = ['[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(160.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #128',
    isomers = [
        '[CH2]C(F)(F)C1(F)[CH]C(F)=C1(590)',
    ],
    reactants = [
        ('CH2CF2(57)', 'FC1=CC(F)=C1(307)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #128',
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

