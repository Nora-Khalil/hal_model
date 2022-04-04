species(
    label = 'F[CH]CC1C[C]1F(9906)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u1 p0 c0 {1,S} {3,S} {4,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (49.7936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,212,367,445,1450,334,575,1197,1424,3202,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27949,0.0518493,-3.21496e-05,7.21841e-09,1.50869e-13,6093.31,25.3517], Tmin=(100,'K'), Tmax=(1222.13,'K')), NASAPolynomial(coeffs=[12.9293,0.023763,-1.00044e-05,1.86258e-09,-1.2907e-13,2495.75,-36.2576], Tmin=(1222.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.7936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'CH2CHF(55)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-153.05,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0219,'amu')),
        NonlinearRotor(inertia=([7.59478,47.6085,55.2033],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([482.184,740.114,880.476,949.569,983.363,1189.2,1343.65,1421.6,1725.69,3171.53,3191.61,3269.97],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09164,-0.0073724,7.45741e-05,-1.12982e-07,5.61696e-11,-18407.3,6.78145], Tmin=(10,'K'), Tmax=(619.705,'K')), NASAPolynomial(coeffs=[1.44203,0.0189088,-1.12569e-05,3.25441e-09,-3.64262e-13,-18255.1,16.8744], Tmin=(619.705,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-153.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CC1(5962)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,D} {7,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (101.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,323,467,575,827,1418,1195.53,1197.63,1198.13,2278.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.17,'J/mol'), sigma=(4.94763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=463.62 K, Pc=55.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08997,-0.0083026,0.000111185,-1.92685e-07,1.08455e-10,12179.8,8.10644], Tmin=(10,'K'), Tmax=(559.269,'K')), NASAPolynomial(coeffs=[1.56308,0.0243173,-1.53208e-05,4.6222e-09,-5.34453e-13,12234.9,16.7948], Tmin=(559.269,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(101.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]1[CH]C1(7380)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (273.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,212,367,445,1450,1756.33,1756.66,1757.17,1757.52,1757.63],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2879.02,'J/mol'), sigma=(5.07614,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=449.70 K, Pc=49.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51151,0.00823174,1.42819e-05,-1.69952e-08,5.05698e-12,32904.9,15.0919], Tmin=(100,'K'), Tmax=(1211.99,'K')), NASAPolynomial(coeffs=[4.17684,0.0157606,-7.07148e-06,1.37064e-09,-9.69304e-14,32029.4,8.80755], Tmin=(1211.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsCsF1s) + radical(cyclopropane)"""),
)

species(
    label = 'F[CH]CC1(F)[CH]C1(11545)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {4,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (39.0531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33005,0.0564073,-4.35415e-05,1.72116e-08,-2.79504e-12,4794.74,24.0026], Tmin=(100,'K'), Tmax=(1425.04,'K')), NASAPolynomial(coeffs=[12.2545,0.0257429,-1.12639e-05,2.11126e-09,-1.45928e-13,1681.2,-32.571], Tmin=(1425.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.0531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs-Cs(F)(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsHH)(H)_ring) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = '[CH2]C(F)C1C[C]1F(9922)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (37.5689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,259,529,569,1128,1321,1390,3140,212,367,445,1450,3000,3100,440,815,1455,1000,180,867.411,867.468,867.475,867.554,867.562,867.573,867.687],'cm^-1')),
        HinderedRotor(inertia=(0.00513298,'amu*angstrom^2'), symmetry=1, barrier=(2.74137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119252,'amu*angstrom^2'), symmetry=1, barrier=(2.74183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31401,0.0503723,-2.50463e-05,-1.52555e-09,3.44688e-12,4622.53,24.1294], Tmin=(100,'K'), Tmax=(1107.38,'K')), NASAPolynomial(coeffs=[12.8174,0.024071,-1.00771e-05,1.89875e-09,-1.33732e-13,1139.75,-36.7629], Tmin=(1107.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH]F(837)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.000263365,8.89188e-06,-1.03032e-08,3.5081e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164132,-7.73117e-07,1.90988e-10,-1.59926e-14,25413.4,-0.815515], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C1C[C]1F(11548)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u1 p0 c0 {2,S} {9,S} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (258.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,212,367,445,1450,3000,3100,440,815,1455,1000,180,1102.07,1102.15,1102.2,1102.33,1102.39,1102.49],'cm^-1')),
        HinderedRotor(inertia=(0.00581937,'amu*angstrom^2'), symmetry=1, barrier=(5.01699,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52044,0.0254094,7.7046e-06,-2.55146e-08,1.09858e-11,31179.4,17.9358], Tmin=(100,'K'), Tmax=(990.777,'K')), NASAPolynomial(coeffs=[8.57029,0.0178555,-6.40066e-06,1.16275e-09,-8.21758e-14,29152.6,-15.3739], Tmin=(990.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Isobutyl)"""),
)

species(
    label = 'FC1CC2CC12F(11549)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {4,S} {6,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-272.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80887,0.0368179,1.00563e-05,-3.44733e-08,1.3987e-11,-32670.3,15.9641], Tmin=(100,'K'), Tmax=(1059.49,'K')), NASAPolynomial(coeffs=[11.6671,0.0257288,-1.12375e-05,2.20314e-09,-1.59925e-13,-36225.8,-39.0866], Tmin=(1059.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FCCC1=C(F)C1(11550)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {2,S} {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-161.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44652,0.0469242,-1.5962e-05,-1.11029e-08,6.9622e-12,-19270.8,23.7518], Tmin=(100,'K'), Tmax=(1051.66,'K')), NASAPolynomial(coeffs=[12.4445,0.0237983,-9.65612e-06,1.81191e-09,-1.2828e-13,-22618.4,-34.7787], Tmin=(1051.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cs-Cd(C)-Cd(F))"""),
)

species(
    label = 'FC=CC1CC1F(9899)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {7,D} {12,S}
7  C u0 p0 c0 {2,S} {6,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-246.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79419,0.0368587,1.15667e-05,-3.82437e-08,1.61205e-11,-29512.1,21.6822], Tmin=(100,'K'), Tmax=(1018.59,'K')), NASAPolynomial(coeffs=[11.8838,0.0243231,-9.86073e-06,1.88665e-09,-1.36394e-13,-32972.7,-34.0776], Tmin=(1018.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cs)"""),
)

species(
    label = 'FCCC1C=C1F(9894)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-140.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60306,0.0545573,-4.3002e-05,1.80606e-08,-3.2211e-12,-16780.1,18.9837], Tmin=(100,'K'), Tmax=(1273.76,'K')), NASAPolynomial(coeffs=[9.60805,0.0294188,-1.33981e-05,2.56625e-09,-1.79986e-13,-18819.3,-21.5727], Tmin=(1273.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = '[CH2]C(F)=CC[CH]F(8699)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u1 p0 c0 {2,S} {3,S} {11,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-110.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,271,519,563,612,1379,334,575,1197,1424,3202,3000,3100,440,815,1455,1000,286.262,1421.67],'cm^-1')),
        HinderedRotor(inertia=(0.0960623,'amu*angstrom^2'), symmetry=1, barrier=(5.57071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327527,'amu*angstrom^2'), symmetry=1, barrier=(19.3192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.654141,'amu*angstrom^2'), symmetry=1, barrier=(38.3291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22956,0.0609087,-5.3773e-05,2.49852e-08,-4.78619e-12,-13156.4,24.7084], Tmin=(100,'K'), Tmax=(1226.31,'K')), NASAPolynomial(coeffs=[11.6885,0.0267934,-1.20436e-05,2.29952e-09,-1.61374e-13,-15721.6,-27.8835], Tmin=(1226.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCsCdF) + radical(Csj(Cs-CdHH)(F1s)(H)) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = '[CH2][CH]F(839)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (116.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0021411,'amu*angstrom^2'), symmetry=1, barrier=(6.24533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56892,0.0121147,-4.46308e-06,3.53001e-10,6.36801e-14,13988.6,11.2355], Tmin=(100,'K'), Tmax=(2106.76,'K')), NASAPolynomial(coeffs=[7.95459,0.00674552,-2.7461e-06,4.76057e-10,-2.99989e-14,11484.3,-14.7486], Tmin=(2106.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(Cs_P)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,25474.2,-0.444973], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]CC1=C(F)C1(11551)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {1,S} {3,S} {5,D}
7  C u1 p0 c0 {2,S} {4,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (32.6752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,334,575,1197,1424,3202,180,180,813.838,813.871,813.886,813.904,813.912],'cm^-1')),
        HinderedRotor(inertia=(0.00775455,'amu*angstrom^2'), symmetry=1, barrier=(3.64503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00776095,'amu*angstrom^2'), symmetry=1, barrier=(3.64737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49572,0.0468345,-2.30243e-05,-2.83291e-09,4.04826e-12,4027.11,24.849], Tmin=(100,'K'), Tmax=(1071.18,'K')), NASAPolynomial(coeffs=[12.4477,0.0212797,-8.72305e-06,1.63727e-09,-1.15617e-13,800.624,-32.8493], Tmin=(1071.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.6752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCsFHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cs-Cd(C)-Cd(F)) + radical(Csj(Cs-CdHH)(F1s)(H))"""),
)

species(
    label = 'F[CH]CC1C=C1F(9902)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u1 p0 c0 {2,S} {4,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (51.8858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,334,575,1197,1424,3202,375.836,380.67,381.501,381.816,2128.34,2129.83,2143.41],'cm^-1')),
        HinderedRotor(inertia=(0.0758966,'amu*angstrom^2'), symmetry=1, barrier=(7.93764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282935,'amu*angstrom^2'), symmetry=1, barrier=(28.4586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51542,0.0592135,-6.5539e-05,4.29126e-08,-1.20059e-11,6325.91,20.4559], Tmin=(100,'K'), Tmax=(848.735,'K')), NASAPolynomial(coeffs=[7.64117,0.0303438,-1.4517e-05,2.83617e-09,-2.01235e-13,5286.07,-8.09264], Tmin=(848.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.8858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(CsCsFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cd(F)-Cd) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'FC=CC1C[C]1F(11552)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {6,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-55.6483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,212,367,445,1450,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,631.611,1293.78,1293.8,1293.89,1294.1,1294.19,1294.27],'cm^-1')),
        HinderedRotor(inertia=(0.231387,'amu*angstrom^2'), symmetry=1, barrier=(5.32005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.09,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82042,0.0409862,-1.63505e-05,-3.09634e-09,2.66445e-12,-6608.79,24.3], Tmin=(100,'K'), Tmax=(1196.56,'K')), NASAPolynomial(coeffs=[10.3733,0.0243423,-1.04633e-05,1.96829e-09,-1.37183e-13,-9510.9,-22.0715], Tmin=(1196.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.6483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)-Cs-Cs) + radical(CsCsCsF1s)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[CH]C[C]1CC1F(8713)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
5  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (2.10893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78901,0.0531686,-5.53308e-05,4.51327e-08,-1.74468e-11,328.995,23.9293], Tmin=(100,'K'), Tmax=(706.453,'K')), NASAPolynomial(coeffs=[3.55469,0.0396118,-1.89879e-05,3.70441e-09,-2.62152e-13,168.344,16.6532], Tmin=(706.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.10893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Tertalkyl) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'F[CH]CC1[CH]C1F(9893)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {4,S} {12,S}
7  C u1 p0 c0 {2,S} {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (55.1133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,180,180,180,806.79,806.86,806.87,806.896,807.061],'cm^-1')),
        HinderedRotor(inertia=(0.00937934,'amu*angstrom^2'), symmetry=1, barrier=(4.33254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00937888,'amu*angstrom^2'), symmetry=1, barrier=(4.33015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3145.12,'J/mol'), sigma=(5.5051,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=491.26 K, Pc=42.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09432,0.0552513,-4.03138e-05,1.45663e-08,-2.10327e-12,6740.41,26.0283], Tmin=(100,'K'), Tmax=(1631.45,'K')), NASAPolynomial(coeffs=[15.2858,0.0204562,-8.32205e-06,1.49326e-09,-9.9971e-14,2109.89,-49.3838], Tmin=(1631.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.1133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'FC[CH]C1C[C]1F(11553)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  C u1 p0 c0 {3,S} {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (52.2296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56155,0.0456804,-1.89666e-05,-3.17435e-09,2.9607e-12,6376.07,25.7565], Tmin=(100,'K'), Tmax=(1197.7,'K')), NASAPolynomial(coeffs=[11.6149,0.025845,-1.13332e-05,2.15538e-09,-1.51194e-13,2982.38,-28.6731], Tmin=(1197.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Cs_S)"""),
)

species(
    label = 'F[CH][CH]C1CC1F(9908)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {7,S} {12,S}
7  C u1 p0 c0 {2,S} {6,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (14.6293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7529,0.0475317,-2.82125e-05,7.68437e-09,-8.25305e-13,1841.46,25.1133], Tmin=(100,'K'), Tmax=(2099.8,'K')), NASAPolynomial(coeffs=[15.6464,0.0210654,-9.30614e-06,1.68178e-09,-1.10642e-13,-3993.26,-52.2217], Tmin=(2099.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.6293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)-Cs-Cs) + radical(Cs_S) + radical(Csj(Cs-CsHH)(F1s)(H))"""),
)

species(
    label = 'FCC[C]1C[C]1F(11554)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {4,S} {7,S}
7  C u1 p0 c0 {2,S} {4,S} {6,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (39.7092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10676,0.0455554,-2.6137e-05,7.0478e-09,-7.73938e-13,4839.39,22.7181], Tmin=(100,'K'), Tmax=(1921.22,'K')), NASAPolynomial(coeffs=[10.5715,0.0279317,-1.23772e-05,2.27313e-09,-1.52629e-13,1586.88,-23.6463], Tmin=(1921.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.7092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring)"""),
)

species(
    label = 'FCCC1[CH][C]1F(9910)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {3,S} {7,S}
7  C u1 p0 c0 {3,S} {6,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (92.7135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,528,1116,1182,1331,1402,1494,3075,3110,212,367,445,1450,529.606,529.608,529.608,529.609,529.61,529.612,529.615,529.617],'cm^-1')),
        HinderedRotor(inertia=(0.000601025,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00060102,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.098,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30622,0.0489006,-1.64217e-05,-1.38612e-08,8.59034e-12,11256.8,25.2059], Tmin=(100,'K'), Tmax=(1032.11,'K')), NASAPolynomial(coeffs=[13.8696,0.0219605,-8.87888e-06,1.68453e-09,-1.20836e-13,7505.01,-41.4142], Tmin=(1032.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.7135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
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
    E0 = (-22.9868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (136.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (222.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (400.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-14.7025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (40.4133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (40.4133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1.98643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (12.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (151.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (173.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (190.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (61.6569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (83.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (316.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (324.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (151.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (143.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (133.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (118.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (132.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (106.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['CH2CHF(55)', 'FC1=CC1(5962)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['F[CH]CC1(F)[CH]C1(11545)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['[CH2]C(F)C1C[C]1F(9922)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]F(837)', '[CH2]C1C[C]1F(11548)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['FC1CC2CC12F(11549)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['FCCC1=C(F)C1(11550)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['FC=CC1CC1F(9899)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['FCCC1C=C1F(9894)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)=CC[CH]F(8699)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]F(839)', 'FC1=CC1(5962)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(6.72436,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'F[CH]CC1=C(F)C1(11551)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(1.92795,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'F[CH]CC1C=C1F(9902)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1020.43,'m^3/(mol*s)'), n=1.51234, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08797594653728344, var=6.974501284715509, Tref=1000.0, N=21, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CHF(55)', 'F[C]1[CH]C1(7380)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(14.0628,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'FC=CC1C[C]1F(11552)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.42983,'m^3/(mol*s)'), n=1.92039, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.021512250841106063, var=1.4639131859069563, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]F(839)', 'F[C]1[CH]C1(7380)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHF(40)', '[CH2]C1C[C]1F(11548)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['F[CH]C[C]1CC1F(8713)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.93363e+09,'s^-1'), n=1.033, Ea=(174.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_Cs2] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]CC1[CH]C1F(9893)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['FC[CH]C1C[C]1F(11553)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(110021,'s^-1'), n=2.37375, Ea=(156.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] + [R2H_S;C_rad_out_1H;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['F[CH][CH]C1CC1F(9908)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.90025e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[CH]CC1C[C]1F(9906)'],
    products = ['FCC[C]1C[C]1F(11554)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(91.367,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FCCC1[CH][C]1F(9910)'],
    products = ['F[CH]CC1C[C]1F(9906)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8619.36,'s^-1'), n=2.2814, Ea=(86.525,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4175',
    isomers = [
        'F[CH]CC1C[C]1F(9906)',
    ],
    reactants = [
        ('CH2CHF(55)', 'FC1=CC1(5962)'),
        ('CH2CHF(55)', 'F[C]1[CH]C1(7380)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4175',
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

