species(
    label = '[CH]=C(C)C([CH2])C(1151)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u1 p0 c0 {1,S} {14,S} {15,S}
6  C u1 p0 c0 {4,D} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (369.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0913524,'amu*angstrom^2'), symmetry=1, barrier=(5.22807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0935813,'amu*angstrom^2'), symmetry=1, barrier=(5.22319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0931744,'amu*angstrom^2'), symmetry=1, barrier=(5.22731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390616,'amu*angstrom^2'), symmetry=1, barrier=(21.6434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01976,0.0586178,-3.51251e-05,6.72555e-09,1.23172e-12,44497.7,25.7158], Tmin=(100,'K'), Tmax=(1046.42,'K')), NASAPolynomial(coeffs=[11.4005,0.0304178,-1.11593e-05,1.94212e-09,-1.30434e-13,41696.6,-27.8401], Tmin=(1046.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C#CC(211)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (172.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,743.639],'cm^-1')),
        HinderedRotor(inertia=(0.630823,'amu*angstrom^2'), symmetry=1, barrier=(14.5039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30524,0.0109264,1.31992e-05,-2.25168e-08,8.87439e-12,20738.2,7.1916], Tmin=(100,'K'), Tmax=(969.981,'K')), NASAPolynomial(coeffs=[5.80052,0.0113536,-4.03494e-06,7.1907e-10,-5.02191e-14,19750,-7.36949], Tmin=(969.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC(60)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31911,0.00817973,3.34731e-05,-4.36187e-08,1.5821e-11,749.326,9.54029], Tmin=(100,'K'), Tmax=(983.759,'K')), NASAPolynomial(coeffs=[5.3676,0.0170742,-6.35103e-06,1.16618e-09,-8.27612e-14,-487.157,-4.54493], Tmin=(983.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(S)(24)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56195e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.5454e-10,-1.62957e-14,50691.8,6.78382], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C(C)C[CH2](1119)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u1 p0 c0 {3,D} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (396.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.135255,'amu*angstrom^2'), symmetry=1, barrier=(3.10979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135252,'amu*angstrom^2'), symmetry=1, barrier=(3.10971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574996,'amu*angstrom^2'), symmetry=1, barrier=(13.2203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3074.1,'J/mol'), sigma=(5.55822,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.17 K, Pc=40.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66998,0.0455229,-2.93065e-05,9.70563e-09,-1.3183e-12,47808.4,21.3998], Tmin=(100,'K'), Tmax=(1683.63,'K')), NASAPolynomial(coeffs=[11.4826,0.0222096,-8.53563e-06,1.4809e-09,-9.70035e-14,44504.3,-31.0522], Tmin=(1683.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(370)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u1 p0 c0 {3,D} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (403.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(3.49973,'amu*angstrom^2'), symmetry=1, barrier=(80.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693833,'amu*angstrom^2'), symmetry=1, barrier=(15.9526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00226253,'amu*angstrom^2'), symmetry=1, barrier=(15.9471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3076.04,'J/mol'), sigma=(5.56648,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.47 K, Pc=40.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91873,0.0379712,-2.3736e-06,-2.30871e-08,1.17455e-11,48660.5,21.7322], Tmin=(100,'K'), Tmax=(953.444,'K')), NASAPolynomial(coeffs=[10.1216,0.0224179,-7.57613e-06,1.29766e-09,-8.83976e-14,46239,-21.9462], Tmin=(953.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(C)=CCC(1240)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {1,S} {4,D} {15,S}
6  C u2 p0 c0 {4,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (279.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0272,0.0571387,-2.35742e-05,-7.50334e-10,2.11832e-12,33679.6,24.4326], Tmin=(100,'K'), Tmax=(1238.1,'K')), NASAPolynomial(coeffs=[10.5796,0.0383695,-1.54849e-05,2.78249e-09,-1.87918e-13,30387.4,-27.4352], Tmin=(1238.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C)C[CH]C(1149)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {1,S} {3,S} {15,S}
6  C u1 p0 c0 {4,D} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (362.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69405,0.0540654,-3.34381e-05,1.15e-08,-1.80014e-12,43640.7,24.4885], Tmin=(100,'K'), Tmax=(1323.25,'K')), NASAPolynomial(coeffs=[6.5627,0.0393482,-1.67551e-05,3.095e-09,-2.12187e-13,42352.2,-0.36354], Tmin=(1323.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = 'CH2(T)(17)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(C)=CC(1186)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  C u2 p0 c0 {3,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (301.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,320.169,322.461,323.839],'cm^-1')),
        HinderedRotor(inertia=(0.691097,'amu*angstrom^2'), symmetry=1, barrier=(51.97,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.688431,'amu*angstrom^2'), symmetry=1, barrier=(51.9658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.687277,'amu*angstrom^2'), symmetry=1, barrier=(52.0417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70885,0.0447872,-1.90929e-05,1.73626e-09,5.4024e-13,36376.5,19.0805], Tmin=(100,'K'), Tmax=(1449.68,'K')), NASAPolynomial(coeffs=[9.4095,0.0309362,-1.24148e-05,2.18487e-09,-1.44101e-13,33366.5,-23.611], Tmin=(1449.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CC1=CCC1C(1241)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (73.6604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87856,0.0290286,5.63868e-05,-9.03667e-08,3.56468e-11,8951.51,18.9776], Tmin=(100,'K'), Tmax=(970.025,'K')), NASAPolynomial(coeffs=[12.7503,0.0280653,-9.95827e-06,1.85063e-09,-1.3526e-13,4778.48,-43.7796], Tmin=(970.025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.6604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C(C)C(=C)C(1242)',
    structure = adjacencyList("""1  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,D} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {13,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (18.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96923,0.0572293,-2.47171e-05,-7.31536e-09,6.63366e-12,2371.22,19.382], Tmin=(100,'K'), Tmax=(1020.73,'K')), NASAPolynomial(coeffs=[13.1747,0.0285538,-1.07259e-05,1.9311e-09,-1.33747e-13,-1118.34,-44.6406], Tmin=(1020.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH3(18)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=C(C)C=C(1121)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {5,D}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  C u1 p0 c0 {2,D} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (303.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.896938,'amu*angstrom^2'), symmetry=1, barrier=(20.6224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894579,'amu*angstrom^2'), symmetry=1, barrier=(20.5681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76568,0.0404502,-9.85315e-06,-1.81201e-08,1.0526e-11,36582.1,17.0182], Tmin=(100,'K'), Tmax=(967.444,'K')), NASAPolynomial(coeffs=[12.2564,0.0175176,-5.992e-06,1.06038e-09,-7.45696e-14,33595.6,-38.1902], Tmin=(967.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(212)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (490.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0859242,'amu*angstrom^2'), symmetry=1, barrier=(1.97557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22084,0.0161542,-7.04072e-06,1.20991e-09,-3.01638e-14,59071.4,10.7688], Tmin=(100,'K'), Tmax=(1807.72,'K')), NASAPolynomial(coeffs=[6.96329,0.0100049,-3.70699e-06,6.32786e-10,-4.05637e-14,57370,-10.4656], Tmin=(1807.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'H(6)',
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
    label = '[CH]=C(C)C(=C)C(1243)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  C u1 p0 c0 {4,D} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
"""),
    E0 = (265.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.739471,'amu*angstrom^2'), symmetry=1, barrier=(17.0019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739504,'amu*angstrom^2'), symmetry=1, barrier=(17.0027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.738945,'amu*angstrom^2'), symmetry=1, barrier=(16.9898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1356,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870805,0.0610352,-4.23219e-05,1.12185e-08,2.6582e-13,32092,20.8319], Tmin=(100,'K'), Tmax=(1057.92,'K')), NASAPolynomial(coeffs=[13.5797,0.0254479,-9.53734e-06,1.69631e-09,-1.15922e-13,28705.5,-44.4935], Tmin=(1057.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C(195)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,S} {7,S}
3 C u1 p0 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.2686e-06,-1.44876e-08,5.7495e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93907,0.0171892,-6.69154e-06,1.21547e-09,-8.39799e-14,33151.2,-4.14876], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C#CC([CH2])C(376)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,T}
5  C u0 p0 c0 {4,T} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (321.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.4595,'amu*angstrom^2'), symmetry=1, barrier=(10.5648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.61196,'amu*angstrom^2'), symmetry=1, barrier=(83.046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0670891,'amu*angstrom^2'), symmetry=1, barrier=(82.9601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05268,0.0371707,-7.10583e-06,-1.96903e-08,1.19936e-11,38774.1,18.6599], Tmin=(100,'K'), Tmax=(877.396,'K')), NASAPolynomial(coeffs=[9.62982,0.0193968,-5.38945e-06,7.89684e-10,-4.88611e-14,36799,-20.5829], Tmin=(877.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(C)=C(C)C(1244)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
2  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {6,S}
6  C u2 p0 c0 {5,S} {16,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (262.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,325,375,415,465,420,450,1700,1750,378.69,378.745,378.786],'cm^-1')),
        HinderedRotor(inertia=(0.528749,'amu*angstrom^2'), symmetry=1, barrier=(53.8235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528656,'amu*angstrom^2'), symmetry=1, barrier=(53.8232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528721,'amu*angstrom^2'), symmetry=1, barrier=(53.8237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528772,'amu*angstrom^2'), symmetry=1, barrier=(53.8242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35627,0.0594785,-3.2876e-05,8.98232e-09,-1.0271e-12,31685.1,21.0651], Tmin=(100,'K'), Tmax=(1826.53,'K')), NASAPolynomial(coeffs=[10.6886,0.0390411,-1.60922e-05,2.85636e-09,-1.88632e-13,28275.9,-29.58], Tmin=(1826.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)=C([CH2])C(1245)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {6,S}
5  C u1 p0 c0 {3,S} {13,S} {14,S}
6  C u1 p0 c0 {4,S} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (194.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889054,0.0597985,-3.35785e-05,4.73205e-09,1.53686e-12,23570.7,21.5241], Tmin=(100,'K'), Tmax=(1127.35,'K')), NASAPolynomial(coeffs=[12.6999,0.030786,-1.21322e-05,2.19515e-09,-1.50418e-13,20088.3,-40.5063], Tmin=(1127.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])C(1246)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {6,D}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {15,S} {16,S}
6  C u0 p0 c0 {3,D} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
"""),
    E0 = (273.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15643,0.0518604,-7.34924e-06,-2.73017e-08,1.45367e-11,32998.8,24.4092], Tmin=(100,'K'), Tmax=(961.967,'K')), NASAPolynomial(coeffs=[13.0072,0.027945,-9.60453e-06,1.66812e-09,-1.14756e-13,29545.3,-38.4036], Tmin=(961.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)C(1247)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {6,D}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {13,S} {14,S}
6  C u0 p0 c0 {3,D} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (327.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18305,0.0559408,-2.80603e-05,-1.88183e-09,5.10357e-12,39438.2,25.815], Tmin=(100,'K'), Tmax=(943.757,'K')), NASAPolynomial(coeffs=[10.5513,0.0305481,-1.04511e-05,1.74948e-09,-1.15378e-13,37032.5,-22.2162], Tmin=(943.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(=C)C(C)C(1248)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  C u2 p0 c0 {4,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
"""),
    E0 = (287.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10222,0.0529928,-4.2515e-06,-2.56624e-08,1.19343e-11,34696.6,23.6521], Tmin=(100,'K'), Tmax=(1031.5,'K')), NASAPolynomial(coeffs=[11.209,0.0368956,-1.44281e-05,2.62104e-09,-1.81466e-13,31382.9,-31.3763], Tmin=(1031.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (74.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (521.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (528.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (269.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (234.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (388.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (83.1073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (138.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (177.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (227.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (202.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (185.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (199.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (481.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (171.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (266.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (204.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (119.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (157.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['C#CC(211)', 'C=CC(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH]=C(C)C[CH2](1119)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(980214,'m^3/(mol*s)'), n=0.201831, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.043461663341785264, var=0.18543765059589312, Tref=1000.0, N=2, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-7R!H=6R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(24)', '[CH]=CC([CH2])C(370)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.55364e+07,'m^3/(mol*s)'), n=-0.357533, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_Sp-6R!H-4CbCdCt',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_Sp-6R!H-4CbCdCt"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['[CH]C(C)=CCC(1240)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['[CH]=C(C)C[CH]C(1149)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(17)', '[CH]C(C)=CC(1186)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['CC1=CCC1C(1241)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['C=C(C)C(=C)C(1242)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(18)', '[CH]=C(C)C=C(1121)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0132,'m^3/(mol*s)'), n=2.41, Ea=(31.9749,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Sp-2R!H=1R!H_N-2R!H->S_4R!H-u0_N-4R!H-inRing_Ext-4R!H-R_Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Sp-2R!H=1R!H_N-2R!H->S_4R!H-u0_N-4R!H-inRing_Ext-4R!H-R_Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]C(212)', 'C=CC(60)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(161.758,'m^3/(mol*s)'), n=1.13767, Ea=(24.5142,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.011839139166081196, var=0.8671401158115553, Tref=1000.0, N=160, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(6)', '[CH]=C(C)C(=C)C(1243)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1348.96,'m^3/(mol*s)'), n=1.432, Ea=(19.4697,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC(211)', '[CH2][CH]C(195)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00205635,'m^3/(mol*s)'), n=2.60319, Ea=(22.5475,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.39655521203258803, var=1.3615091655633467, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_N-Sp-2R!H=1R!H_Ext-3C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_N-Sp-2R!H=1R!H_Ext-3C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(18)', 'C#CC([CH2])C(376)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.102,'m^3/(mol*s)'), n=2.41, Ea=(36.1485,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing_Ext-4R!H-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing_Ext-4R!H-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]C(212)', '[CH2][CH]C(195)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.65983e+11,'m^3/(mol*s)'), n=-0.84129, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.1366575604485066, var=2.914302663837648, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_N-3BrCO->O_Ext-2CF-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_N-3BrCO->O_Ext-2CF-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C(C)=C(C)C(1244)'],
    products = ['[CH]=C(C)C([CH2])C(1151)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.614e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['[CH2]C(C)=C([CH2])C(1245)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['[CH2]C(=C)C([CH2])C(1246)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['[CH2]C([CH2])C(=C)C(1247)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C)C([CH2])C(1151)'],
    products = ['[CH]C(=C)C(C)C(1248)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #685',
    isomers = [
        '[CH]=C(C)C([CH2])C(1151)',
    ],
    reactants = [
        ('C#CC(211)', 'C=CC(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #685',
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

