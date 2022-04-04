species(
    label = '[O]CC([O])O(7113)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {9,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-165.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,2144.9,2144.95],'cm^-1')),
        HinderedRotor(inertia=(0.164454,'amu*angstrom^2'), symmetry=1, barrier=(3.78113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163787,'amu*angstrom^2'), symmetry=1, barrier=(3.76579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42524,0.0451959,-8.59069e-05,9.4247e-08,-3.78173e-11,-19824.8,21.0516], Tmin=(100,'K'), Tmax=(848.592,'K')), NASAPolynomial(coeffs=[-0.435689,0.0323366,-1.66082e-05,3.22024e-09,-2.22259e-13,-18390.7,39.9732], Tmin=(848.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'CH2O(20)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CO(1075)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-389.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,910.228,910.24],'cm^-1')),
        HinderedRotor(inertia=(0.000203447,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76621,-0.000867468,3.3288e-05,-3.95813e-08,1.3994e-11,-46830.2,7.56895], Tmin=(100,'K'), Tmax=(997.627,'K')), NASAPolynomial(coeffs=[5.52472,0.00779112,-3.35078e-06,6.868e-10,-5.23895e-14,-47962.8,-4.82863], Tmin=(997.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,29230.2,5.12617], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C[CH]O(803)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {8,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-10.3955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,180,2261.05],'cm^-1')),
        HinderedRotor(inertia=(0.223718,'amu*angstrom^2'), symmetry=1, barrier=(5.14371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223266,'amu*angstrom^2'), symmetry=1, barrier=(5.13333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47466,0.0416324,-7.70031e-05,7.78854e-08,-2.9046e-11,-1203.08,16.406], Tmin=(100,'K'), Tmax=(877.572,'K')), NASAPolynomial(coeffs=[2.30007,0.0226909,-1.08908e-05,2.03334e-09,-1.36524e-13,-412.424,21.5557], Tmin=(877.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.3955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])O(7019)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {8,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-19.5078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1935.18],'cm^-1')),
        HinderedRotor(inertia=(0.2007,'amu*angstrom^2'), symmetry=1, barrier=(4.61448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20003,'amu*angstrom^2'), symmetry=1, barrier=(4.59908,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66061,0.0330653,-4.66659e-05,4.29179e-08,-1.60973e-11,-2301.52,17.2788], Tmin=(100,'K'), Tmax=(804.621,'K')), NASAPolynomial(coeffs=[3.98157,0.0198211,-9.52741e-06,1.83298e-09,-1.27458e-13,-2297.94,12.5363], Tmin=(804.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.5078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = 'OC1COO1(7117)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {4,S} {9,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-215.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61752,0.0293608,-1.72906e-05,4.78393e-09,-5.28722e-13,-25921.6,13.2054], Tmin=(100,'K'), Tmax=(2016.23,'K')), NASAPolynomial(coeffs=[10.1472,0.014423,-6.17775e-06,1.10956e-09,-7.31337e-14,-28958,-28.4013], Tmin=(2016.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxetane)"""),
)

species(
    label = 'O=C(O)CO(2273)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-588.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34774,0.0307996,-1.42115e-05,-2.66903e-09,2.84857e-12,-70757.7,17.268], Tmin=(100,'K'), Tmax=(1080.1,'K')), NASAPolynomial(coeffs=[9.50826,0.014602,-6.04928e-06,1.13943e-09,-8.05089e-14,-72906.5,-20.6158], Tmin=(1080.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-588.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs)"""),
)

species(
    label = 'O=CC(O)O(7186)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-537.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42658,0.0336613,-2.77733e-05,1.1974e-08,-2.12183e-12,-64626.4,18.8422], Tmin=(100,'K'), Tmax=(1319.95,'K')), NASAPolynomial(coeffs=[8.61062,0.0149213,-6.47722e-06,1.21813e-09,-8.46773e-14,-66258.9,-12.7089], Tmin=(1319.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]CC=O(802)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-90.1149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,584.194],'cm^-1')),
        HinderedRotor(inertia=(0.069063,'amu*angstrom^2'), symmetry=1, barrier=(16.7259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04362,0.0136821,1.69878e-05,-3.2263e-08,1.32993e-11,-10797.3,13.3581], Tmin=(100,'K'), Tmax=(973.607,'K')), NASAPolynomial(coeffs=[8.68703,0.0086086,-3.10042e-06,5.99582e-10,-4.54029e-14,-12754.6,-18.1256], Tmin=(973.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.1149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][O](132)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88412,-0.00363928,3.2856e-05,-4.13634e-08,1.59641e-11,23210.8,7.4797], Tmin=(100,'K'), Tmax=(933.047,'K')), NASAPolynomial(coeffs=[6.69324,0.000290178,8.61305e-07,-1.56324e-10,7.33557e-15,21991.3,-9.60369], Tmin=(933.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
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
    label = '[O]CC(=O)O(2355)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {8,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-345.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180,704.754,704.756,704.776,704.792],'cm^-1')),
        HinderedRotor(inertia=(0.00834502,'amu*angstrom^2'), symmetry=1, barrier=(2.94145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00834161,'amu*angstrom^2'), symmetry=1, barrier=(2.94032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53281,0.0284229,-1.94174e-05,6.58115e-09,-8.99583e-13,-41451.6,17.3916], Tmin=(100,'K'), Tmax=(1695.16,'K')), NASAPolynomial(coeffs=[9.47927,0.0120315,-4.91306e-06,8.76888e-10,-5.83221e-14,-43806.6,-19.7872], Tmin=(1695.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-345.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]O(1357)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (14.9737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2183.26],'cm^-1')),
        HinderedRotor(inertia=(0.0143998,'amu*angstrom^2'), symmetry=1, barrier=(48.7076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36704,0.0314251,-9.8653e-05,1.28233e-07,-5.44488e-11,1806.86,9.28624], Tmin=(100,'K'), Tmax=(862.563,'K')), NASAPolynomial(coeffs=[-6.65905,0.0315132,-1.81054e-05,3.60548e-09,-2.49747e-13,5262.84,66.1811], Tmin=(862.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.9737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: thermo_DFT_CCSDTF12_BAC + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = '[O]C(O)C=O(7164)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {8,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-294.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3550.66],'cm^-1')),
        HinderedRotor(inertia=(0.051076,'amu*angstrom^2'), symmetry=1, barrier=(2.09405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625001,'amu*angstrom^2'), symmetry=1, barrier=(14.37,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73016,0.0299169,-2.84184e-05,1.57732e-08,-3.78631e-12,-35325.3,19.2328], Tmin=(100,'K'), Tmax=(971.494,'K')), NASAPolynomial(coeffs=[6.01113,0.0164081,-7.56105e-06,1.46047e-09,-1.032e-13,-35962.8,3.49885], Tmin=(971.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C[C](O)O(7187)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {9,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-185.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11693,0.0484539,-8.3142e-05,7.99449e-08,-2.93661e-11,-22270.9,21.2939], Tmin=(100,'K'), Tmax=(850.731,'K')), NASAPolynomial(coeffs=[4.06425,0.0235997,-1.164e-05,2.22268e-09,-1.52214e-13,-22034.1,15.5529], Tmin=(850.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C(O)[CH]O(7188)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-210.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,286.814,2138.02],'cm^-1')),
        HinderedRotor(inertia=(0.105969,'amu*angstrom^2'), symmetry=1, barrier=(6.18722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106168,'amu*angstrom^2'), symmetry=1, barrier=(6.18587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106103,'amu*angstrom^2'), symmetry=1, barrier=(6.18363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83645,0.0561396,-0.000102692,9.8248e-08,-3.52743e-11,-25262.8,21.8557], Tmin=(100,'K'), Tmax=(869.154,'K')), NASAPolynomial(coeffs=[4.79615,0.0228584,-1.13245e-05,2.14161e-09,-1.44837e-13,-25034.7,12.2638], Tmin=(869.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C(O)O(7189)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-210.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83645,0.0561396,-0.000102692,9.8248e-08,-3.52743e-11,-25262.8,21.1625], Tmin=(100,'K'), Tmax=(869.154,'K')), NASAPolynomial(coeffs=[4.79615,0.0228584,-1.13245e-05,2.14161e-09,-1.44837e-13,-25034.7,11.5707], Tmin=(869.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][C](O)CO(7190)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-185.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11693,0.0484539,-8.3142e-05,7.99449e-08,-2.93661e-11,-22270.9,21.9871], Tmin=(100,'K'), Tmax=(850.731,'K')), NASAPolynomial(coeffs=[4.06425,0.0235997,-1.164e-05,2.22268e-09,-1.52214e-13,-22034.1,16.246], Tmin=(850.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C([O])CO(7191)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {9,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-165.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,180,2144.82,2145.08],'cm^-1')),
        HinderedRotor(inertia=(0.164204,'amu*angstrom^2'), symmetry=1, barrier=(3.77537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164024,'amu*angstrom^2'), symmetry=1, barrier=(3.77123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42524,0.0451959,-8.59069e-05,9.4247e-08,-3.78173e-11,-19824.8,20.3584], Tmin=(100,'K'), Tmax=(848.592,'K')), NASAPolynomial(coeffs=[-0.435689,0.0323366,-1.66082e-05,3.22024e-09,-2.22259e-13,-18390.7,39.2801], Tmin=(848.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]COO(7193)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {9,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {3,S} {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (61.1095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,180,1219.49],'cm^-1')),
        HinderedRotor(inertia=(0.268047,'amu*angstrom^2'), symmetry=1, barrier=(6.16292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267792,'amu*angstrom^2'), symmetry=1, barrier=(6.15707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.35725,'amu*angstrom^2'), symmetry=1, barrier=(54.1978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77798,0.0588236,-0.000109184,1.07843e-07,-3.98746e-11,7420.23,19.6676], Tmin=(100,'K'), Tmax=(859.737,'K')), NASAPolynomial(coeffs=[3.55203,0.0274102,-1.39694e-05,2.67781e-09,-1.82874e-13,7971.1,16.3548], Tmin=(859.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.1095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH)"""),
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
    E0 = (-87.2158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (310.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (301.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-78.9315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-23.8156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-23.8156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (48.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-42.6666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-17.4364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-13.9099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (11.4596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (285.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (70.8206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-3.33706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-11.9351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-11.9351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (18.4256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (246.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['CH2O(20)', 'O=CO(1075)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[O]C[CH]O(803)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[CH2]C([O])O(7019)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['OC1COO1(7117)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['O=C(O)CO(2273)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['O=CC(O)O(7186)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(7)', '[O]CC=O(802)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(32.072,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][O](132)', 'O=CO(1075)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(75.9091,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]CC(=O)O(2355)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(37.8687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH2O(20)', '[O][CH]O(1357)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(59.332,'m^3/(mol*s)'), n=1.2975, Ea=(12.1674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.2890935075667379, var=2.146502499593098, Tref=1000.0, N=39, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]C(O)C=O(7164)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.99,'m^3/(mol*s)'), n=2.12, Ea=(15.7295,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_N-Sp-5R!H=4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_N-Sp-5R!H=4C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][O](132)', '[O][CH]O(1357)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['[O]C[C](O)O(7187)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(O)[CH]O(7188)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['[O][CH]C(O)O(7189)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]CC([O])O(7113)'],
    products = ['[O][C](O)CO(7190)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C([O])CO(7191)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(188.023,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][CH]COO(7193)'],
    products = ['[O]CC([O])O(7113)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

network(
    label = 'PDepNetwork #2777',
    isomers = [
        '[O]CC([O])O(7113)',
    ],
    reactants = [
        ('CH2O(20)', 'O=CO(1075)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2777',
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

