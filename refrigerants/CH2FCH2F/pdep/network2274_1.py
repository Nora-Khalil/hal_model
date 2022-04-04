species(
    label = 'C=C(CC=O)OF(6002)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {3,D} {4,S} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-149.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,725.021,738.218],'cm^-1')),
        HinderedRotor(inertia=(0.00989997,'amu*angstrom^2'), symmetry=1, barrier=(3.77505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812598,'amu*angstrom^2'), symmetry=1, barrier=(18.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.811162,'amu*angstrom^2'), symmetry=1, barrier=(18.6502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57694,0.0498379,-3.44328e-05,1.06413e-08,-1.28362e-12,-17938.1,22.9064], Tmin=(100,'K'), Tmax=(1917.68,'K')), NASAPolynomial(coeffs=[17.3672,0.0169016,-8.67011e-06,1.68504e-09,-1.16027e-13,-23994.2,-63.5538], Tmin=(1917.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.56331e-09,3.13594e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88611e-07,1.21037e-10,-7.84012e-15,-14180.9,6.71043], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=C(C)OF(2585)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-54.0769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.164664,'amu*angstrom^2'), symmetry=1, barrier=(3.78595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0340763,'amu*angstrom^2'), symmetry=1, barrier=(7.51179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0696,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74756,0.0275218,-9.67407e-06,-2.80673e-09,1.90981e-12,-6504.87,9.64687], Tmin=(10,'K'), Tmax=(1181.53,'K')), NASAPolynomial(coeffs=[7.44792,0.0216801,-1.07455e-05,2.58688e-09,-2.44738e-13,-7845.96,-10.7974], Tmin=(1181.53,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-54.0769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""CDC(C)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CCC(=O)CF(6084)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {3,D} {4,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-491.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73924,0.0544114,-5.62367e-05,3.54027e-08,-9.79896e-12,-59047.6,22.9293], Tmin=(100,'K'), Tmax=(845.812,'K')), NASAPolynomial(coeffs=[6.64593,0.0312066,-1.50838e-05,2.96573e-09,-2.11354e-13,-59877.6,0.0790539], Tmin=(845.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-491.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + group(Cds-OdCsH)"""),
)

species(
    label = 'C=COC(=C)OF(6000)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {4,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {2,S} {7,D} {8,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-29.5168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,391.575,391.575,391.576,391.576],'cm^-1')),
        HinderedRotor(inertia=(0.0871667,'amu*angstrom^2'), symmetry=1, barrier=(9.48433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0871663,'amu*angstrom^2'), symmetry=1, barrier=(9.48432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246088,'amu*angstrom^2'), symmetry=1, barrier=(26.7761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.82,'J/mol'), sigma=(5.32487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.02 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23798,0.0641398,-7.40386e-05,4.70788e-08,-1.22818e-11,-3453.47,23.219], Tmin=(100,'K'), Tmax=(923.309,'K')), NASAPolynomial(coeffs=[10.0452,0.0259851,-1.20534e-05,2.32339e-09,-1.63719e-13,-5079.84,-18.5681], Tmin=(923.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2sCF) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'FO[C]1C[CH]OC1(6085)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {1,S} {6,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u1 p0 c0 {2,S} {4,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (124.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.77561,0.0578537,-5.39409e-05,2.70574e-08,-5.1673e-12,15045.4,22.1753], Tmin=(100,'K'), Tmax=(1505.42,'K')), NASAPolynomial(coeffs=[12.8287,0.0162152,-2.87432e-06,2.01256e-10,-3.01174e-15,12505.7,-37.2864], Tmin=(1505.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + radical(C2CsJO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(C=C[O])OF(6086)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  C u0 p0 c0 {3,S} {5,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (89.4384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,271.463,271.468,271.469],'cm^-1')),
        HinderedRotor(inertia=(0.398315,'amu*angstrom^2'), symmetry=1, barrier=(20.8293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398315,'amu*angstrom^2'), symmetry=1, barrier=(20.8293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398303,'amu*angstrom^2'), symmetry=1, barrier=(20.8293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660423,0.0616505,-4.11877e-05,-2.82793e-09,8.65422e-12,10887.9,27.4058], Tmin=(100,'K'), Tmax=(960.179,'K')), NASAPolynomial(coeffs=[19.8927,0.00925549,-2.64707e-06,4.84427e-10,-3.79176e-14,5916.58,-71.2523], Tmin=(960.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.4384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = 'C=C([CH]C[O])OF(6087)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
5  C u1 p0 c0 {4,S} {6,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (109.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,317.093,317.3,1468.27,1469.61],'cm^-1')),
        HinderedRotor(inertia=(0.148729,'amu*angstrom^2'), symmetry=1, barrier=(10.6229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618102,'amu*angstrom^2'), symmetry=1, barrier=(44.2522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622599,'amu*angstrom^2'), symmetry=1, barrier=(44.2392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49517,0.063582,-6.40369e-05,1.47646e-08,1.69718e-11,13254.1,22.4444], Tmin=(100,'K'), Tmax=(531.808,'K')), NASAPolynomial(coeffs=[7.35802,0.0330959,-1.64402e-05,3.22519e-09,-2.27886e-13,12438,-3.94789], Tmin=(531.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C(C[C]=O)OF(6088)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u1 p0 c0 {3,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (106.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950,354.583,355.191],'cm^-1')),
        HinderedRotor(inertia=(0.0838426,'amu*angstrom^2'), symmetry=1, barrier=(7.42025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0838701,'amu*angstrom^2'), symmetry=1, barrier=(7.46303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0831918,'amu*angstrom^2'), symmetry=1, barrier=(7.41315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155587,'amu*angstrom^2'), symmetry=1, barrier=(13.8506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752778,0.0765606,-0.000114437,9.31076e-08,-3.01981e-11,12919.1,27.2747], Tmin=(100,'K'), Tmax=(811.639,'K')), NASAPolynomial(coeffs=[10.2109,0.0255382,-1.19921e-05,2.267e-09,-1.55632e-13,11529,-15.4867], Tmin=(811.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'C[C](C[C]=O)OF(6089)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (74.9993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1855,455,950,203.844,2205.58],'cm^-1')),
        HinderedRotor(inertia=(0.299373,'amu*angstrom^2'), symmetry=1, barrier=(8.82725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299371,'amu*angstrom^2'), symmetry=1, barrier=(8.8273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299378,'amu*angstrom^2'), symmetry=1, barrier=(8.82729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299351,'amu*angstrom^2'), symmetry=1, barrier=(8.82731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883404,0.0777173,-0.000130424,1.1891e-07,-4.17469e-11,9123.75,26.5905], Tmin=(100,'K'), Tmax=(856.056,'K')), NASAPolynomial(coeffs=[6.7481,0.0310964,-1.50608e-05,2.84549e-09,-1.93484e-13,8823.81,3.32102], Tmin=(856.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.9993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C(CC[O])OF(6090)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u1 p0 c0 {6,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (239.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3120,650,792.5,1650,180,1258.64,1260.03],'cm^-1')),
        HinderedRotor(inertia=(0.268126,'amu*angstrom^2'), symmetry=1, barrier=(6.16475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270108,'amu*angstrom^2'), symmetry=1, barrier=(6.2103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268726,'amu*angstrom^2'), symmetry=1, barrier=(6.17853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12341,0.0708716,-0.000110412,1.00209e-07,-3.61065e-11,28925.4,26.893], Tmin=(100,'K'), Tmax=(819.287,'K')), NASAPolynomial(coeffs=[5.98379,0.0333388,-1.64232e-05,3.16523e-09,-2.1935e-13,28592.3,7.24045], Tmin=(819.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C[CH]O)OF(6091)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,S} {11,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {4,S} {10,S}
7  C u1 p0 c0 {5,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (194.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,363.366,363.383],'cm^-1')),
        HinderedRotor(inertia=(0.00127648,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108669,'amu*angstrom^2'), symmetry=1, barrier=(10.1853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100777,'amu*angstrom^2'), symmetry=1, barrier=(9.44246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100796,'amu*angstrom^2'), symmetry=1, barrier=(9.44235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.529667,0.0818702,-0.000127359,1.04342e-07,-3.3564e-11,23487.6,27.715], Tmin=(100,'K'), Tmax=(838.253,'K')), NASAPolynomial(coeffs=[11.2592,0.0237841,-1.10943e-05,2.07575e-09,-1.41016e-13,21930.8,-20.7123], Tmin=(838.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(C=CO)OF(6092)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {6,S} {12,S}
4  C u0 p0 c0 {2,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {8,S}
6  C u0 p0 c0 {3,S} {5,D} {9,S}
7  C u0 p0 c0 {4,D} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-170.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0559,'amu*angstrom^2'), symmetry=1, barrier=(24.2773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05815,'amu*angstrom^2'), symmetry=1, barrier=(24.3291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05568,'amu*angstrom^2'), symmetry=1, barrier=(24.2721,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377496,0.0771566,-8.53645e-05,4.38586e-08,-8.30564e-12,-20387.5,25.2809], Tmin=(100,'K'), Tmax=(1500.29,'K')), NASAPolynomial(coeffs=[22.4185,0.00341571,1.32298e-06,-4.20964e-10,3.25539e-14,-25768.7,-89.0815], Tmin=(1500.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
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
    label = 'C=C([O])CC=O(6093)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {2,D} {3,S} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-169.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,244.852,4000],'cm^-1')),
        HinderedRotor(inertia=(0.561672,'amu*angstrom^2'), symmetry=1, barrier=(23.8707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.560811,'amu*angstrom^2'), symmetry=1, barrier=(23.8706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0812,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91842,0.0361418,-4.85413e-06,-1.84164e-08,8.73512e-12,-20350.2,20.7004], Tmin=(100,'K'), Tmax=(1074.82,'K')), NASAPolynomial(coeffs=[12.4652,0.0177008,-8.15958e-06,1.64714e-09,-1.21496e-13,-23819.4,-36.5345], Tmin=(1074.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]F(127)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=[C]CC=O(6094)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,D} {2,S} {8,S}
4  C u0 p0 c0 {5,D} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {4,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (144.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,594.479],'cm^-1')),
        HinderedRotor(inertia=(0.0326171,'amu*angstrom^2'), symmetry=1, barrier=(16.3861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706247,'amu*angstrom^2'), symmetry=1, barrier=(16.238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0818,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5693,0.0252417,3.86036e-06,-1.67085e-08,6.0059e-12,17460.3,17.2901], Tmin=(100,'K'), Tmax=(1209.11,'K')), NASAPolynomial(coeffs=[8.86207,0.0205062,-1.02165e-05,2.05383e-09,-1.4805e-13,14763,-19.1252], Tmin=(1209.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OF(196)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (235.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2950,3100,1380,975,1025,1650,1685,370,1858.56],'cm^-1')),
        HinderedRotor(inertia=(0.190472,'amu*angstrom^2'), symmetry=1, barrier=(8.45471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27376,0.0120477,-4.21846e-06,6.03973e-11,9.31488e-14,28267.1,10.4556], Tmin=(100,'K'), Tmax=(2667.03,'K')), NASAPolynomial(coeffs=[19.9694,-0.00396407,5.52783e-07,-7.38518e-11,6.52156e-15,17217.5,-85.6831], Tmin=(2667.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=CJO)"""),
)

species(
    label = 'CH2CHO(35)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,530.062,530.609,975.335,1640.89,1642.3],'cm^-1')),
        HinderedRotor(inertia=(0.00115918,'amu*angstrom^2'), symmetry=1, barrier=(2.21299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C(=C)OF(2610)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {8,S} {9,S}
5 C u0 p0 c0 {3,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (109.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,343.549],'cm^-1')),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(12.8688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0184306,'amu*angstrom^2'), symmetry=1, barrier=(33.8385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0616,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87336,0.00851378,0.000122333,-3.04869e-07,2.24407e-10,13187.1,10.7552], Tmin=(10,'K'), Tmax=(466.39,'K')), NASAPolynomial(coeffs=[4.30812,0.0285537,-1.8564e-05,5.80221e-09,-6.95957e-13,12888,6.21793], Tmin=(466.39,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(109.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[CH2]C(DC)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(C=C[O])OF(6095)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {8,S}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {5,D} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-29.4956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17495,'amu*angstrom^2'), symmetry=1, barrier=(27.0143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17453,'amu*angstrom^2'), symmetry=1, barrier=(27.0049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.290357,0.0671264,-7.10495e-05,3.55575e-08,-6.67486e-12,-3401.57,24.2964], Tmin=(100,'K'), Tmax=(1460.25,'K')), NASAPolynomial(coeffs=[19.8141,0.00623537,-8.88726e-07,5.07802e-11,-9.99556e-16,-8313.4,-74.5809], Tmin=(1460.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.4956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C=C(C[C]=O)OF(6096)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (10.0694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,452.627,453.077],'cm^-1')),
        HinderedRotor(inertia=(0.0008209,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101397,'amu*angstrom^2'), symmetry=1, barrier=(14.8148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102511,'amu*angstrom^2'), symmetry=1, barrier=(14.823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76852,0.0500234,-4.12652e-05,1.62438e-08,-2.58803e-12,1290.13,23.0144], Tmin=(100,'K'), Tmax=(1447.69,'K')), NASAPolynomial(coeffs=[12.4582,0.0204877,-1.06624e-05,2.15113e-09,-1.54391e-13,-1804.94,-32.5119], Tmin=(1447.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.0694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C(CC=O)OF(6097)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {3,D} {4,S} {10,S}
7  C u1 p0 c0 {5,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (97.205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,536.035],'cm^-1')),
        HinderedRotor(inertia=(0.0727748,'amu*angstrom^2'), symmetry=1, barrier=(1.67324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0896777,'amu*angstrom^2'), symmetry=1, barrier=(18.2412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0888939,'amu*angstrom^2'), symmetry=1, barrier=(18.2411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87023,0.0494216,-3.89773e-05,1.43911e-08,-2.15575e-12,11764.9,22.2294], Tmin=(100,'K'), Tmax=(1514,'K')), NASAPolynomial(coeffs=[12.3072,0.021847,-1.16574e-05,2.36121e-09,-1.69284e-13,8604.58,-32.4514], Tmin=(1514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C#COF(428)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (221.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([374,480,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.63561,'amu*angstrom^2'), symmetry=1, barrier=(37.6059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9094,0.00583843,6.25691e-05,-1.62148e-07,1.16725e-10,26695.5,7.9678], Tmin=(10,'K'), Tmax=(505.498,'K')), NASAPolynomial(coeffs=[5.5535,0.011202,-7.86683e-06,2.64783e-09,-3.37635e-13,26294.6,-1.16432], Tmin=(505.498,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(221.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C#COF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=CO(461)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {7,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u0 p0 c0 {2,D} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H3OH""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]C(CC=O)OF(6098)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {3,D} {4,S} {11,S}
7  C u0 p1 c0 {5,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (186.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.914516,0.0754977,-0.000119432,1.07137e-07,-3.82049e-11,22478.4,25.3797], Tmin=(100,'K'), Tmax=(809.205,'K')), NASAPolynomial(coeffs=[7.14731,0.0322634,-1.62579e-05,3.16268e-09,-2.20352e-13,21876.5,-0.856864], Tmin=(809.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsCsOsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-147.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (138.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-26.1638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (35.3544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (138.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (70.7515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (90.8333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (128.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (67.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (223.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (177.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (5.85852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-125.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (205.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (194.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (102.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (146.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (182.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (267.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (116.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (156.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(CC=O)OF(6002)'],
    products = ['HF(38)', 'CH2CO(28)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(44.2884,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'C=C(C)OF(2585)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.53038,'m^3/(mol*s)'), n=2.0364, Ea=(352.59,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.26314699591391894, var=4.003648838577825, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(CC=O)OF(6002)'],
    products = ['O=CCC(=O)CF(6084)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(165.276,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=COC(=C)OF(6000)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(106.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['FO[C]1C[CH]OC1(6085)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.43612e+09,'s^-1'), n=0.0758668, Ea=(56.4791,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(C=C[O])OF(6086)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([CH]C[O])OF(6087)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(C[C]=O)OF(6088)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[C](C[C]=O)OF(6089)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(CC[O])OF(6090)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C[CH]O)OF(6091)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(C=CO)OF(6092)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(75.1006,'s^-1'), n=3.11041, Ea=(218.365,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00848789959541425, var=6.520606768146176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R!H->C',), comment="""Estimated from node Root_3R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'C=C([O])CC=O(6093)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(13.5053,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(127)', 'C=[C]CC=O(6094)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OF(196)', 'CH2CHO(35)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HCO(15)', '[CH2]C(=C)OF(2610)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.23652e+08,'m^3/(mol*s)'), n=-0.382134, Ea=(1.83252,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3183737345414065, var=0.1827168304887792, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-2CF-R_Sp-4R!H=3BrBrCCOO',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-2CF-R_Sp-4R!H=3BrBrCCOO
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'C=C(C=C[O])OF(6095)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(5.69885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'C=C(C[C]=O)OF(6096)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(2.19974,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[CH]=C(CC=O)OF(6097)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -5.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(CC=O)OF(6002)'],
    products = ['C#COF(428)', 'C=CO(461)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.56726e+11,'s^-1'), n=0.169307, Ea=(308.231,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_5R!H->O_Ext-2C-R_N-7R!H->C',), comment="""Estimated from node Root_1R!H->C_2R!H->C_5R!H->O_Ext-2C-R_N-7R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(CC=O)OF(6098)'],
    products = ['C=C(CC=O)OF(6002)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(11.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

network(
    label = 'PDepNetwork #2274',
    isomers = [
        'C=C(CC=O)OF(6002)',
    ],
    reactants = [
        ('HF(38)', 'CH2CO(28)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2274',
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

