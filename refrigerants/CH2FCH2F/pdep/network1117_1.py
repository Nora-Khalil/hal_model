species(
    label = 'C#COC(O)[CH]F(2455)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u1 p0 c0 {1,S} {4,S} {9,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-166.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.953272,'amu*angstrom^2'), symmetry=1, barrier=(21.9176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.954075,'amu*angstrom^2'), symmetry=1, barrier=(21.9361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952545,'amu*angstrom^2'), symmetry=1, barrier=(21.9009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.954191,'amu*angstrom^2'), symmetry=1, barrier=(21.9387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332053,0.0823781,-0.000114976,7.73198e-08,-2.00811e-11,-19942.1,23.2304], Tmin=(100,'K'), Tmax=(951.234,'K')), NASAPolynomial(coeffs=[16.6993,0.013551,-6.44014e-06,1.25102e-09,-8.85071e-14,-23055.9,-54.9135], Tmin=(951.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Csj(Cs-O2sO2sH)(F1s)(H))"""),
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
    label = 'O=C[CH]F(435)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-191.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,1215,1347,1486,3221,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.21522,'amu*angstrom^2'), symmetry=1, barrier=(35.2035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2974.27,'J/mol'), sigma=(4.91649,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.57 K, Pc=56.79 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96722,0.00209998,6.04576e-05,-1.2409e-07,8.01311e-11,-23018.4,8.77531], Tmin=(10,'K'), Tmax=(478.869,'K')), NASAPolynomial(coeffs=[2.63637,0.0192853,-1.23829e-05,3.78104e-09,-4.41625e-13,-22960.5,13.4894], Tmin=(478.869,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-191.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]C(O)OC#C-2(6630)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {3,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
4  C u2 p0 c0 {3,S} {9,S}
5  C u0 p0 c0 {1,S} {6,T}
6  C u0 p0 c0 {5,T} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (271.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2175,525,750,770,3400,2100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.447009,0.0656224,-7.48591e-05,3.92443e-08,-7.67038e-12,32829.1,23.0601], Tmin=(100,'K'), Tmax=(1391.2,'K')), NASAPolynomial(coeffs=[20.143,0.00247073,2.62733e-07,-1.23779e-10,9.58648e-15,27979.9,-76.1965], Tmin=(1391.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(223)',
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
    label = 'C#CO[CH]O(3046)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u1 p0 c0 {1,S} {2,S} {6,S}
4 C u0 p0 c0 {1,S} {5,T}
5 C u0 p0 c0 {4,T} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (58.3708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1311,'amu*angstrom^2'), symmetry=1, barrier=(26.0062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13152,'amu*angstrom^2'), symmetry=1, barrier=(26.016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12872,'amu*angstrom^2'), symmetry=1, barrier=(25.9514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801234,0.0569919,-7.09583e-05,3.89699e-08,-7.74763e-12,7147.07,17.4889], Tmin=(100,'K'), Tmax=(1448.61,'K')), NASAPolynomial(coeffs=[18.5345,-0.00345231,3.51528e-06,-7.73383e-10,5.5179e-14,3213.68,-70.4789], Tmin=(1448.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.3708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(OCJO)"""),
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
    label = 'C#COC(O)[C]F(3849)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u2 p0 c0 {1,S} {4,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (5.32689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2175,525,750,770,3400,2100,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.508321,0.0724503,-9.0896e-05,5.2955e-08,-1.17429e-11,770.209,22.3883], Tmin=(100,'K'), Tmax=(1118.92,'K')), NASAPolynomial(coeffs=[18.8486,0.00688717,-3.00491e-06,5.89292e-10,-4.30264e-14,-3334.12,-68.154], Tmin=(1118.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.32689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CsCCl_triplet)"""),
)

species(
    label = 'OC1O[C]=CC1F(6631)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {4,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-289.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88719,0.0348249,8.28654e-06,-4.25725e-08,2.11192e-11,-34688.8,18.5079], Tmin=(100,'K'), Tmax=(923.892,'K')), NASAPolynomial(coeffs=[14.267,0.0112998,-2.34534e-06,3.3179e-10,-2.41545e-14,-38259.8,-47.184], Tmin=(923.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1OC(O)C1F(6632)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u1 p0 c0 {6,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-208.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22366,0.0481028,-1.30862e-05,-3.1661e-08,1.99606e-11,-24965,18.1273], Tmin=(100,'K'), Tmax=(912.923,'K')), NASAPolynomial(coeffs=[18.871,0.00528726,5.65809e-07,-2.27118e-10,1.4459e-14,-29625.1,-73.2786], Tmin=(912.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsCCFH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cds_P)"""),
)

species(
    label = 'HCCO(21)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (166.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87608,0.0221205,-3.58869e-05,3.05403e-08,-1.01281e-11,20163.4,13.6968], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91479,0.00371409,-1.30137e-06,2.06473e-10,-1.21477e-14,19359.6,-5.50567], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(166.705,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'OC=CF(455)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-317.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(1.37732,'amu*angstrom^2'), symmetry=1, barrier=(31.6673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95033,0.00268027,6.73468e-05,-1.15374e-07,5.90231e-11,-38204.1,8.04413], Tmin=(10,'K'), Tmax=(637.71,'K')), NASAPolynomial(coeffs=[1.89968,0.0278733,-2.09144e-05,7.21532e-09,-9.21467e-13,-38193.3,15.049], Tmin=(637.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-317.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OH(7)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#COC=CF(2007)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {5,S}
3 C u0 p0 c0 {2,S} {4,D} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {8,S}
5 C u0 p0 c0 {2,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (14.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48797,'amu*angstrom^2'), symmetry=1, barrier=(34.2113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48927,'amu*angstrom^2'), symmetry=1, barrier=(34.2412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3150.18,'J/mol'), sigma=(5.12424,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.05 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87804,0.0456264,-4.5193e-05,2.2646e-08,-4.54094e-12,1778.28,17.5753], Tmin=(100,'K'), Tmax=(1198.5,'K')), NASAPolynomial(coeffs=[11.0529,0.0150052,-6.86843e-06,1.32774e-09,-9.40678e-14,-420.913,-28.3491], Tmin=(1198.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'C#COC(O)=CF(6633)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {1,S} {4,D} {8,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-142.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,194,682,905,1196,1383,3221,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.829705,'amu*angstrom^2'), symmetry=1, barrier=(19.0766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832976,'amu*angstrom^2'), symmetry=1, barrier=(19.1517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45978,'amu*angstrom^2'), symmetry=1, barrier=(33.5632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.952804,0.0717471,-0.000109959,8.74238e-08,-2.75117e-11,-17092.9,20.6862], Tmin=(100,'K'), Tmax=(780.796,'K')), NASAPolynomial(coeffs=[11.1365,0.0195862,-9.77122e-06,1.89662e-09,-1.32305e-13,-18683.4,-25.9263], Tmin=(780.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'O[CH][CH]F(456)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u1 p0 c0 {2,S} {4,S} {5,S}
4 C u1 p0 c0 {1,S} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-62.0342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.78444,'amu*angstrom^2'), symmetry=1, barrier=(18.0358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783638,'amu*angstrom^2'), symmetry=1, barrier=(18.0174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06762,0.0405792,-5.09762e-05,3.02879e-08,-6.85208e-12,-7389.71,15.1943], Tmin=(100,'K'), Tmax=(1097.19,'K')), NASAPolynomial(coeffs=[11.964,0.00450027,-1.65186e-06,3.17801e-10,-2.32669e-14,-9561.35,-33.4679], Tmin=(1097.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.0342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJOH) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = 'C2H(22)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C(O)[CH]F(2543)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {1,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-207.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.892082,'amu*angstrom^2'), symmetry=1, barrier=(20.5107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0695885,'amu*angstrom^2'), symmetry=1, barrier=(41.2584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93343,0.0512034,-8.54769e-05,8.01951e-08,-2.99492e-11,-24885.5,17.2154], Tmin=(100,'K'), Tmax=(773.954,'K')), NASAPolynomial(coeffs=[5.94451,0.0218506,-1.18773e-05,2.4036e-09,-1.7155e-13,-25248.2,0.560417], Tmin=(773.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-O2sO2sH)(F1s)(H))"""),
)

species(
    label = 'C#CO[CH][CH]F(2038)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {5,S}
3 C u1 p0 c0 {2,S} {4,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {8,S}
5 C u0 p0 c0 {2,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (209.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20668,'amu*angstrom^2'), symmetry=1, barrier=(27.7439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20013,'amu*angstrom^2'), symmetry=1, barrier=(27.5934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20539,'amu*angstrom^2'), symmetry=1, barrier=(27.7144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625165,0.0699392,-9.1363e-05,5.24505e-08,-1.05769e-11,25303.3,19.6068], Tmin=(100,'K'), Tmax=(924.79,'K')), NASAPolynomial(coeffs=[18.9396,0.00215798,3.09161e-08,-6.41221e-11,5.20302e-15,21426.9,-69.9615], Tmin=(924.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOCs) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = 'C#COC([O])[CH]F(6634)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u1 p0 c0 {1,S} {4,S} {9,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (58.8104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,334,575,1197,1424,3202,2175,525,750,770,3400,2100,238.22,238.261,238.385,238.404],'cm^-1')),
        HinderedRotor(inertia=(0.583979,'amu*angstrom^2'), symmetry=1, barrier=(23.5174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583888,'amu*angstrom^2'), symmetry=1, barrier=(23.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583376,'amu*angstrom^2'), symmetry=1, barrier=(23.5168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587608,0.0809776,-0.000132472,1.10065e-07,-3.59747e-11,7190.51,23.1682], Tmin=(100,'K'), Tmax=(778.135,'K')), NASAPolynomial(coeffs=[12.0535,0.0198918,-1.05819e-05,2.09241e-09,-1.46829e-13,5471.06,-28.8546], Tmin=(778.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.8104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(Csj(Cs-O2sO2sH)(F1s)(H))"""),
)

species(
    label = 'C#CO[C](O)[CH]F(6635)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u1 p0 c0 {2,S} {3,S} {5,S}
5  C u1 p0 c0 {1,S} {4,S} {8,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (38.3516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180.002,181.078,181.135],'cm^-1')),
        HinderedRotor(inertia=(0.894464,'amu*angstrom^2'), symmetry=1, barrier=(20.8168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905614,'amu*angstrom^2'), symmetry=1, barrier=(20.8219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905552,'amu*angstrom^2'), symmetry=1, barrier=(20.8205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905549,'amu*angstrom^2'), symmetry=1, barrier=(20.8204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.339564,0.0834627,-0.000126664,9.12036e-08,-2.5235e-11,4741.86,23.8913], Tmin=(100,'K'), Tmax=(895.349,'K')), NASAPolynomial(coeffs=[16.5058,0.0112384,-5.66285e-06,1.10661e-09,-7.77687e-14,1847,-52.3144], Tmin=(895.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.3516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cs_P) + radical(Csj(Cs-O2sO2sH)(F1s)(H))"""),
)

species(
    label = '[C]#COC(O)[CH]F(6636)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u1 p0 c0 {1,S} {4,S} {9,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u1 p0 c0 {6,T}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (170.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,334,575,1197,1424,3202,2175,525,237.931,238.043,238.045,238.281],'cm^-1')),
        HinderedRotor(inertia=(0.478032,'amu*angstrom^2'), symmetry=1, barrier=(19.2285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478268,'amu*angstrom^2'), symmetry=1, barrier=(19.2253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478052,'amu*angstrom^2'), symmetry=1, barrier=(19.2262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.477877,'amu*angstrom^2'), symmetry=1, barrier=(19.2276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199339,0.0878527,-0.000144771,1.14264e-07,-3.44894e-11,20609.1,24.4423], Tmin=(100,'K'), Tmax=(876.641,'K')), NASAPolynomial(coeffs=[15.4986,0.012306,-5.68639e-06,1.02668e-09,-6.71261e-14,18147.1,-46.096], Tmin=(876.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Csj(Cs-O2sO2sH)(F1s)(H)) + radical(Acetyl)"""),
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
    label = 'C#COC(O)[C]F-2(3858)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u0 p1 c0 {1,S} {4,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (1.24111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,617,898,1187,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.021,'amu*angstrom^2'), symmetry=1, barrier=(23.4748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02076,'amu*angstrom^2'), symmetry=1, barrier=(23.4692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02068,'amu*angstrom^2'), symmetry=1, barrier=(23.4674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02072,'amu*angstrom^2'), symmetry=1, barrier=(23.4683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.588221,0.0685117,-8.05486e-05,4.25527e-08,-8.06396e-12,277.882,23.7932], Tmin=(100,'K'), Tmax=(1009.63,'K')), NASAPolynomial(coeffs=[18.9786,0.00559275,-1.83964e-06,3.32785e-10,-2.4476e-14,-3942.28,-67.6152], Tmin=(1009.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.24111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CJ2_singlet-FCs) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CO[C](O)CF(1392)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {5,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u1 p0 c0 {2,S} {3,S} {4,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-160.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,551,1088,1226,1380,1420,1481,3057,3119,360,370,350,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05606,'amu*angstrom^2'), symmetry=1, barrier=(24.2809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05689,'amu*angstrom^2'), symmetry=1, barrier=(24.3,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05664,'amu*angstrom^2'), symmetry=1, barrier=(24.2941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05672,'amu*angstrom^2'), symmetry=1, barrier=(24.2961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345235,0.0730818,-8.5601e-05,4.70787e-08,-9.8612e-12,-19105.8,24.0452], Tmin=(100,'K'), Tmax=(1184.65,'K')), NASAPolynomial(coeffs=[19.6941,0.00775039,-2.87928e-06,5.27094e-10,-3.7371e-14,-23690.1,-72.5807], Tmin=(1184.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cs_P)"""),
)

species(
    label = 'C#COC([O])CF(1405)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u0 p0 c0 {6,T} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-139.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,2175,525,750,770,3400,2100,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2029,'amu*angstrom^2'), symmetry=1, barrier=(27.657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2036,'amu*angstrom^2'), symmetry=1, barrier=(27.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20313,'amu*angstrom^2'), symmetry=1, barrier=(27.6623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03211,0.0652798,-7.2052e-05,3.94841e-08,-8.5687e-12,-16676.1,21.7566], Tmin=(100,'K'), Tmax=(1117.36,'K')), NASAPolynomial(coeffs=[14.1926,0.0181674,-8.80644e-06,1.74941e-09,-1.25944e-13,-19617.1,-43.1958], Tmin=(1117.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ)"""),
)

species(
    label = '[C]#COC(O)CF(6637)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {4,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {7,T}
7  C u1 p0 c0 {6,T}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
"""),
    E0 = (-28.1034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,528,1116,1182,1331,1402,1494,3075,3110,2175,525,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.987923,'amu*angstrom^2'), symmetry=1, barrier=(22.7143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.988057,'amu*angstrom^2'), symmetry=1, barrier=(22.7174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.988054,'amu*angstrom^2'), symmetry=1, barrier=(22.7173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.988155,'amu*angstrom^2'), symmetry=1, barrier=(22.7196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448873,0.0745239,-9.30845e-05,5.5953e-08,-1.28924e-11,-3249.04,23.7265], Tmin=(100,'K'), Tmax=(1075.57,'K')), NASAPolynomial(coeffs=[17.6395,0.0105919,-3.92332e-06,6.87814e-10,-4.66656e-14,-6946.93,-60.4603], Tmin=(1075.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.1034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
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
    E0 = (-102.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (294.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (222.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (166.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-158.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-156.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-191.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-7.87331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (18.4001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (54.2595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (299.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (187.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (220.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (199.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (331.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (146.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (167.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-59.7725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-105.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-33.8812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#COC(O)[CH]F(2455)'],
    products = ['CH2CO(28)', 'O=C[CH]F(435)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(114.99,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F(37)', '[CH]C(O)OC#C-2(6630)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]F(223)', 'C#CO[CH]O(3046)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/O2;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(5)', 'C#COC(O)[C]F(3849)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#COC(O)[CH]F(2455)'],
    products = ['OC1O[C]=CC1F(6631)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C#COC(O)[CH]F(2455)'],
    products = ['[CH]=C1OC(O)C1F(6632)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(60.6686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCCO(21)', 'OC=CF(455)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(10.2449,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(7)', 'C#COC=CF(2007)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(565.393,'m^3/(mol*s)'), n=1.28408, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.06277471317878429, var=1.0052381843076323, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C#COC(O)=CF(6633)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(21820,'m^3/(mol*s)'), n=0.859, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['HCCO(21)', 'O[CH][CH]F(456)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C2H(22)', '[O]C(O)[CH]F(2543)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OH(7)', 'C#CO[CH][CH]F(2038)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'C#COC([O])[CH]F(6634)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'C#CO[C](O)[CH]F(6635)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.43711e+23,'m^3/(mol*s)'), n=-5.99271, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25884854535692575, var=10.966526674850428, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl
Ea raised from -1.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[C]#COC(O)[CH]F(6636)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.43711e+23,'m^3/(mol*s)'), n=-5.99271, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25884854535692575, var=10.966526674850428, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHF(40)', 'C#CO[CH]O(3046)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'C#COC(O)[C]F-2(3858)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(5.20252,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CO[C](O)CF(1392)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.81094e+06,'s^-1'), n=1.92799, Ea=(150.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_1H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#COC([O])CF(1405)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(18000,'s^-1'), n=2.287, Ea=(84.6842,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[C]#COC(O)CF(6637)'],
    products = ['C#COC(O)[CH]F(2455)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4708.23,'s^-1'), n=2.0907, Ea=(44.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_1H] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_1H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1117',
    isomers = [
        'C#COC(O)[CH]F(2455)',
    ],
    reactants = [
        ('CH2CO(28)', 'O=C[CH]F(435)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1117',
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

