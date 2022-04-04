species(
    label = 'C=[C]CO[C](C)F(7010)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (42.8435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0907357,'amu*angstrom^2'), symmetry=1, barrier=(2.08619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933915,'amu*angstrom^2'), symmetry=1, barrier=(2.14725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0923063,'amu*angstrom^2'), symmetry=1, barrier=(2.1223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702298,'amu*angstrom^2'), symmetry=1, barrier=(16.1472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03281,0.0711971,-9.39906e-05,7.93806e-08,-2.83935e-11,5254.06,26.498], Tmin=(100,'K'), Tmax=(754.622,'K')), NASAPolynomial(coeffs=[6.21754,0.0384928,-1.86028e-05,3.6097e-09,-2.53339e-13,4620.24,3.92955], Tmin=(754.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.8435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)F(499)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-453.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,486,617,768,1157,1926],'cm^-1')),
        HinderedRotor(inertia=(0.163766,'amu*angstrom^2'), symmetry=1, barrier=(3.76529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3044.49,'J/mol'), sigma=(4.92747,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.54 K, Pc=57.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95746,0.00861016,1.67316e-05,-2.39571e-08,8.75423e-12,-54563.9,8.39893], Tmin=(10,'K'), Tmax=(963.489,'K')), NASAPolynomial(coeffs=[3.63904,0.0176469,-9.3479e-06,2.39871e-09,-2.40789e-13,-54860.6,8.06499], Tmin=(963.489,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-453.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=C=C(6556)',
    structure = adjacencyList("""1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37446,0.00704632,2.78302e-05,-3.99438e-08,1.55726e-11,21188.6,7.62051], Tmin=(100,'K'), Tmax=(949.709,'K')), NASAPolynomial(coeffs=[6.79959,0.00959973,-3.02065e-06,5.37819e-10,-3.92599e-14,19772.3,-12.7584], Tmin=(949.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[C]F-2(838)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u2 p0 c0 {1,S} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (86.2806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1986.18],'cm^-1')),
        HinderedRotor(inertia=(0.0255602,'amu*angstrom^2'), symmetry=1, barrier=(9.94152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51663,0.00840623,6.86137e-06,-1.03526e-08,3.33162e-12,10396.3,8.5806], Tmin=(100,'K'), Tmax=(1175.76,'K')), NASAPolynomial(coeffs=[4.70435,0.0108934,-4.6396e-06,8.90531e-10,-6.30071e-14,9665.83,0.739437], Tmin=(1175.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.2806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = 'C=[C]C[O](689)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {4,D} {7,S} {8,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (319.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.110137,'amu*angstrom^2'), symmetry=1, barrier=(2.53227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0632,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63077,0.0408503,-8.2829e-05,9.16991e-08,-3.63161e-11,38513.3,13.7006], Tmin=(100,'K'), Tmax=(865.187,'K')), NASAPolynomial(coeffs=[-0.468585,0.0281972,-1.41119e-05,2.70319e-09,-1.84571e-13,40059.5,34.0405], Tmin=(865.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[C]=C(134)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u2 p0 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.36689e-08,8.45455e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.85899e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[C](C)F(981)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u1 p0 c0 {2,S} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-91.4372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,3000,3100,440,815,1455,1000,297.905,298.072],'cm^-1')),
        HinderedRotor(inertia=(0.0849761,'amu*angstrom^2'), symmetry=1, barrier=(5.35507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0850337,'amu*angstrom^2'), symmetry=1, barrier=(5.35627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0850769,'amu*angstrom^2'), symmetry=1, barrier=(5.35631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0696,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10969,0.0451222,-5.11282e-05,3.4213e-08,-9.68701e-12,-10932.4,17.8924], Tmin=(100,'K'), Tmax=(842.734,'K')), NASAPolynomial(coeffs=[6.92076,0.0222853,-1.04778e-05,2.05351e-09,-1.46225e-13,-11743.2,-4.49472], Tmin=(842.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.4372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(CsJOCH3)"""),
)

species(
    label = 'C=C1COC1(C)F(7012)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-277.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48214,0.0426611,7.39876e-06,-3.62296e-08,1.53381e-11,-33262.2,18.5678], Tmin=(100,'K'), Tmax=(1046.36,'K')), NASAPolynomial(coeffs=[13.1334,0.0270164,-1.15963e-05,2.26425e-09,-1.64509e-13,-37282.3,-45.7292], Tmin=(1046.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C=COC(C)F(6521)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {12,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-270.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.49669,0.0613391,-1.98407e-05,-2.94484e-08,1.85738e-11,-32409.8,22.2877], Tmin=(100,'K'), Tmax=(958.334,'K')), NASAPolynomial(coeffs=[21.0805,0.0133002,-3.93378e-06,7.27122e-10,-5.66123e-14,-38094.3,-85.2158], Tmin=(958.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CCOC(=C)F(7013)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {7,D}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-261.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894373,0.0594936,-4.10503e-05,1.08062e-08,-3.26095e-14,-31290.3,24.3484], Tmin=(100,'K'), Tmax=(1135.55,'K')), NASAPolynomial(coeffs=[14.5663,0.0235081,-9.5968e-06,1.78133e-09,-1.2423e-13,-35180.3,-46.8046], Tmin=(1135.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][C]=C(6610)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {3,S} {6,S} {7,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (395.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.110742,'amu*angstrom^2'), symmetry=1, barrier=(23.6192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28577,0.0103677,1.67968e-05,-2.70052e-08,1.05232e-11,47575.3,10.4038], Tmin=(100,'K'), Tmax=(980.89,'K')), NASAPolynomial(coeffs=[6.52797,0.0104288,-3.60851e-06,6.6848e-10,-4.85522e-14,46300.3,-8.4326], Tmin=(980.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C[C]([O])F(1344)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-67.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,4000],'cm^-1')),
        HinderedRotor(inertia=(0.378092,'amu*angstrom^2'), symmetry=1, barrier=(8.69308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91286,0.0277959,-4.27461e-05,4.09524e-08,-1.51725e-11,-8085.5,14.0056], Tmin=(100,'K'), Tmax=(854.905,'K')), NASAPolynomial(coeffs=[3.47688,0.0163951,-7.36921e-06,1.37672e-09,-9.36531e-14,-7861.75,13.2456], Tmin=(854.905,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = 'C=C=CO[C](C)F(6531)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {2,S} {7,D} {11,S}
6  C u0 p0 c0 {7,D} {12,S} {13,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-76.0833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07998,'amu*angstrom^2'), symmetry=1, barrier=(24.8309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08294,'amu*angstrom^2'), symmetry=1, barrier=(24.899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07698,'amu*angstrom^2'), symmetry=1, barrier=(24.7619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406028,0.0693462,-5.83767e-05,1.62918e-08,1.30769e-12,-9012.66,23.4159], Tmin=(100,'K'), Tmax=(984.899,'K')), NASAPolynomial(coeffs=[19.0009,0.0140161,-4.85768e-06,8.78881e-10,-6.31344e-14,-13654.7,-70.9817], Tmin=(984.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.0833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'C=[C]COC(=C)F(7016)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  C u0 p0 c0 {7,D} {12,S} {13,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-23.3119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,326,540,652,719,1357,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180,180,705.285],'cm^-1')),
        HinderedRotor(inertia=(0.70687,'amu*angstrom^2'), symmetry=1, barrier=(16.2523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131914,'amu*angstrom^2'), symmetry=1, barrier=(16.2522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706812,'amu*angstrom^2'), symmetry=1, barrier=(16.251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16916,0.0605422,-5.3762e-05,2.48579e-08,-4.68097e-12,-2700.51,23.8309], Tmin=(100,'K'), Tmax=(1257.9,'K')), NASAPolynomial(coeffs=[12.5382,0.0243903,-1.06527e-05,2.01085e-09,-1.40319e-13,-5560.75,-33.6265], Tmin=(1257.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.3119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C#CCO[C](C)F(7017)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {6,T} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-29.0483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,2175,525,750,770,3400,2100,180,180,2755.04],'cm^-1')),
        HinderedRotor(inertia=(0.66473,'amu*angstrom^2'), symmetry=1, barrier=(15.2834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.45652,'amu*angstrom^2'), symmetry=1, barrier=(56.4802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664679,'amu*angstrom^2'), symmetry=1, barrier=(15.2823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4563,'amu*angstrom^2'), symmetry=1, barrier=(56.4752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01334,0.0697103,-8.84054e-05,6.60188e-08,-2.03787e-11,-3389.89,23.7157], Tmin=(100,'K'), Tmax=(785.531,'K')), NASAPolynomial(coeffs=[8.70158,0.0305552,-1.36261e-05,2.54548e-09,-1.74914e-13,-4597.58,-11.5185], Tmin=(785.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.0483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs))"""),
)

species(
    label = 'C[C]F(124)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (80.9091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1049.2],'cm^-1')),
        HinderedRotor(inertia=(0.0264144,'amu*angstrom^2'), symmetry=1, barrier=(3.87782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1971.36,'J/mol'), sigma=(5.118e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50526,0.00692148,1.56045e-05,-2.26115e-08,8.49098e-12,9752.45,9.13904], Tmin=(100,'K'), Tmax=(981.304,'K')), NASAPolynomial(coeffs=[5.2739,0.00958149,-3.54759e-06,6.48822e-10,-4.59773e-14,8930.14,-1.78146], Tmin=(981.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[CH2]C=CO[C](C)F(6534)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {6,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {6,D} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {5,D} {12,S}
7  C u1 p0 c0 {5,S} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-101.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.999844,'amu*angstrom^2'), symmetry=1, barrier=(22.9884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999145,'amu*angstrom^2'), symmetry=1, barrier=(22.9723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9966,'amu*angstrom^2'), symmetry=1, barrier=(22.9138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998136,'amu*angstrom^2'), symmetry=1, barrier=(22.9491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374358,0.0671146,-4.04429e-05,-5.94659e-09,9.83959e-12,-12028.2,24.7565], Tmin=(100,'K'), Tmax=(965.244,'K')), NASAPolynomial(coeffs=[20.0357,0.0149129,-4.81551e-06,8.68549e-10,-6.38914e-14,-17187.6,-76.4674], Tmin=(965.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CsCsF1sO2s) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(F)OC[C]=C(7018)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  C u0 p0 c0 {7,D} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (61.0188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10203,0.0693705,-8.13996e-05,5.92753e-08,-1.86418e-11,7438.1,26.8464], Tmin=(100,'K'), Tmax=(757.839,'K')), NASAPolynomial(coeffs=[7.16822,0.0373515,-1.8023e-05,3.52242e-09,-2.49463e-13,6518.67,-0.737351], Tmin=(757.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.0188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(Cs-F1sO2sH)(H)(H)) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCO[C](C)F(7019)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  C u0 p0 c0 {3,S} {7,D} {13,S}
7  C u1 p0 c0 {6,D} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (52.0978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,395,473,707,1436,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,214.362,214.85,215.66],'cm^-1')),
        HinderedRotor(inertia=(0.00363179,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193421,'amu*angstrom^2'), symmetry=1, barrier=(6.4548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194931,'amu*angstrom^2'), symmetry=1, barrier=(6.45611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610506,'amu*angstrom^2'), symmetry=1, barrier=(20.3458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16779,0.0677828,-7.62809e-05,5.26761e-08,-1.57426e-11,6363,25.9594], Tmin=(100,'K'), Tmax=(794.097,'K')), NASAPolynomial(coeffs=[7.3303,0.0367409,-1.76443e-05,3.44854e-09,-2.44532e-13,5384.29,-2.3503], Tmin=(794.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.0978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=COC(C)F(6538)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {12,S}
6  C u1 p0 c0 {7,S} {13,S} {14,S}
7  C u1 p0 c0 {5,D} {6,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-57.9008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.423154,0.0638134,-2.70744e-05,-2.19102e-08,1.59158e-11,-6821.61,24.2522], Tmin=(100,'K'), Tmax=(959.175,'K')), NASAPolynomial(coeffs=[20.9857,0.013596,-4.11114e-06,7.521e-10,-5.7563e-14,-12400.8,-82.6138], Tmin=(959.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.9008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C](F)OCC=C(7020)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u1 p0 c0 {1,S} {2,S} {7,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (17.7019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36553,0.0660627,-6.29954e-05,1.76723e-08,1.17971e-11,2216.08,26.6041], Tmin=(100,'K'), Tmax=(547.252,'K')), NASAPolynomial(coeffs=[6.92343,0.0375691,-1.81446e-05,3.53908e-09,-2.50043e-13,1426.12,1.48133], Tmin=(547.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.7019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
)

species(
    label = '[CH]=[C]COC(C)F(7021)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
6  C u1 p0 c0 {5,S} {7,D}
7  C u1 p0 c0 {6,D} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (95.4148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,1685,370,3120,650,792.5,1650,291.461,291.461,291.461],'cm^-1')),
        HinderedRotor(inertia=(0.19426,'amu*angstrom^2'), symmetry=1, barrier=(11.7104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00453578,'amu*angstrom^2'), symmetry=1, barrier=(11.7104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194258,'amu*angstrom^2'), symmetry=1, barrier=(11.7104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.505821,'amu*angstrom^2'), symmetry=1, barrier=(30.4917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25801,0.0657236,-6.81685e-05,4.24966e-08,-1.15471e-11,11569.8,24.997], Tmin=(100,'K'), Tmax=(863.863,'K')), NASAPolynomial(coeffs=[7.55205,0.0365798,-1.75632e-05,3.4429e-09,-2.4496e-13,10482.4,-4.44697], Tmin=(863.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.4148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (-24.6517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (338.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (441.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-16.3674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (38.7484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (14.5733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-6.1979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (50.6145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (80.0958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (127.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (119.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (260.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (333.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (114.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (159.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (90.0395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (116.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (117.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (60.9598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['CC(=O)F(499)', 'C=C=C(6556)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[C]F-2(838)', 'C=[C]C[O](689)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C(134)', '[CH2]O[C](C)F(981)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['C=C1COC1(C)F(7012)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['C=C=COC(C)F(6521)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['C=CCOC(=C)F(7013)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(=O)F(499)', '[CH2][C]=C(6610)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(119.641,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[C]([O])F(1344)', 'C=C=C(6556)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.60832e-16,'m^3/(mol*s)'), n=6.08142, Ea=(9.69694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.42494241855262277, var=2.081425725945606, Tref=1000.0, N=380, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C=C=CO[C](C)F(6531)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(11.8696,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'C=[C]COC(=C)F(7016)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.67e+06,'m^3/(mol*s)'), n=0.1, Ea=(6.07236,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_4R!H->O_Sp-4O-2C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_4R!H->O_Sp-4O-2C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C#CCO[C](C)F(7017)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(172591,'m^3/(mol*s)'), n=0.983154, Ea=(3.98016,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[C]([O])F(1344)', '[CH2][C]=C(6610)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[C]F(124)', 'C=[C]C[O](689)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['[CH2]C=CO[C](C)F(6534)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['[CH2]C(F)OC[C]=C(7018)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.465e+11,'s^-1'), n=0, Ea=(184.64,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 440 used for R2H_S;C_rad_out_noH;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CCO[C](C)F(7019)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['[CH2][C]=COC(C)F(6538)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.90025e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]CO[C](C)F(7010)'],
    products = ['[CH2][C](F)OCC=C(7020)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]COC(C)F(7021)'],
    products = ['C=[C]CO[C](C)F(7010)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2401',
    isomers = [
        'C=[C]CO[C](C)F(7010)',
    ],
    reactants = [
        ('CC(=O)F(499)', 'C=C=C(6556)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2401',
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

