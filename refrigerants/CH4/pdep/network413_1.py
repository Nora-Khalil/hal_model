species(
    label = '[CH2]CO[CH]O(784)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4  C u1 p0 c0 {3,S} {8,S} {9,S}
5  C u1 p0 c0 {1,S} {2,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-27.2325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,317.883,317.906,317.914],'cm^-1')),
        HinderedRotor(inertia=(0.00166826,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00166799,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157332,'amu*angstrom^2'), symmetry=1, barrier=(11.2813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157328,'amu*angstrom^2'), symmetry=1, barrier=(11.2814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37404,0.0598781,-7.87403e-05,5.57531e-08,-1.57155e-11,-3182.68,20.9194], Tmin=(100,'K'), Tmax=(869.468,'K')), NASAPolynomial(coeffs=[10.3,0.018814,-7.89691e-06,1.43385e-09,-9.69491e-14,-4734.86,-20.8948], Tmin=(869.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.2325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + radical(CJCO) + radical(OCJO)"""),
)

species(
    label = 'O=CO(77)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76624,-0.000867782,3.32891e-05,-3.95828e-08,1.39946e-11,-46830.2,7.56886], Tmin=(100,'K'), Tmax=(997.616,'K')), NASAPolynomial(coeffs=[5.52461,0.0077913,-3.35088e-06,6.86825e-10,-5.23916e-14,-47962.8,-4.82803], Tmin=(997.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H4(29)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]O(118)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.873,3308.85],'cm^-1')),
        HinderedRotor(inertia=(0.0103145,'amu*angstrom^2'), symmetry=1, barrier=(22.7124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75937,0.00296141,8.9037e-06,-1.3501e-08,5.39791e-12,24775.6,6.76289], Tmin=(100,'K'), Tmax=(940.438,'K')), NASAPolynomial(coeffs=[5.09115,0.00321234,-9.31657e-07,1.59608e-10,-1.15723e-14,24263.5,-0.97116], Tmin=(940.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[O](66)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (188.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1398.33],'cm^-1')),
        HinderedRotor(inertia=(0.00547724,'amu*angstrom^2'), symmetry=1, barrier=(7.58298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57171,0.0102136,5.9092e-06,-7.99876e-09,2.0708e-12,22733,11.7516], Tmin=(100,'K'), Tmax=(1490.84,'K')), NASAPolynomial(coeffs=[4.74094,0.0150201,-6.91918e-06,1.3118e-09,-8.98248e-14,21501.6,2.68328], Tmin=(1490.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCOJ) + radical(CJCO)"""),
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
    label = '[CH2]O[CH]O(64)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u1 p0 c0 {1,S} {2,S} {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-13.4029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,430.052,430.593],'cm^-1')),
        HinderedRotor(inertia=(0.000917461,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0666157,'amu*angstrom^2'), symmetry=1, barrier=(8.66446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0657148,'amu*angstrom^2'), symmetry=1, barrier=(8.64312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19801,0.0396778,-5.03099e-05,3.31155e-08,-8.53807e-12,-1547.15,16.2847], Tmin=(100,'K'), Tmax=(955.291,'K')), NASAPolynomial(coeffs=[9.27221,0.0100563,-3.79755e-06,6.55633e-10,-4.32052e-14,-2898.72,-17.5206], Tmin=(955.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.4029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsHHH) + radical(OCJO) + radical(CsJOCH3)"""),
)

species(
    label = 'OC1CCO1(785)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {11,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-296.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68434,0.0174522,4.5954e-05,-7.73616e-08,3.39722e-11,-35608.5,15.3585], Tmin=(100,'K'), Tmax=(872.741,'K')), NASAPolynomial(coeffs=[10.4133,0.0134602,-1.20734e-06,-6.9531e-11,1.06585e-14,-38154.6,-27.7351], Tmin=(872.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-296.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Oxetane)"""),
)

species(
    label = 'C=COCO(788)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-339.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37662,0.0287175,6.81097e-05,-1.38415e-07,6.39755e-11,-40760.7,18.1825], Tmin=(100,'K'), Tmax=(903.01,'K')), NASAPolynomial(coeffs=[28.5135,-0.0144455,1.18311e-05,-2.38452e-09,1.5767e-13,-48802.9,-127.361], Tmin=(903.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CCOC=O(259)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {2,D} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-407.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49142,0.0261853,8.77816e-06,-2.29355e-08,8.50625e-12,-48965,16.9474], Tmin=(100,'K'), Tmax=(1107.16,'K')), NASAPolynomial(coeffs=[7.6405,0.0239689,-1.04198e-05,1.99229e-09,-1.41067e-13,-51109.5,-12.9537], Tmin=(1107.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsH)"""),
)

species(
    label = '[CH2][CH2](81)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,971.411,1436.01,2688.07,2689.99],'cm^-1')),
        HinderedRotor(inertia=(0.0117579,'amu*angstrom^2'), symmetry=1, barrier=(17.2429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40736,0.0100313,6.40915e-06,-1.41289e-08,5.92664e-12,38288.2,6.11704], Tmin=(100,'K'), Tmax=(954.263,'K')), NASAPolynomial(coeffs=[5.5225,0.00856171,-2.90742e-06,5.02351e-10,-3.4457e-14,37547.7,-5.75282], Tmin=(954.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]COC=O(225)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4  C u1 p0 c0 {3,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-196.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,1285.65],'cm^-1')),
        HinderedRotor(inertia=(0.136405,'amu*angstrom^2'), symmetry=1, barrier=(3.13622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51738,'amu*angstrom^2'), symmetry=1, barrier=(34.8876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0297441,'amu*angstrom^2'), symmetry=1, barrier=(34.8877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31655,0.031598,-1.15536e-05,-3.30113e-09,2.26195e-12,-23511.8,18.2117], Tmin=(100,'K'), Tmax=(1210.43,'K')), NASAPolynomial(coeffs=[8.90154,0.0194875,-8.50484e-06,1.60639e-09,-1.12038e-13,-25812.9,-17.7349], Tmin=(1210.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(CJCO)"""),
)

species(
    label = '[O][CH]O(212)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1989.6],'cm^-1')),
        HinderedRotor(inertia=(0.0936924,'amu*angstrom^2'), symmetry=1, barrier=(2.15417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40205,0.0209569,-4.96195e-05,5.95345e-08,-2.44629e-11,660.882,10.3083], Tmin=(100,'K'), Tmax=(868.742,'K')), NASAPolynomial(coeffs=[-0.343261,0.0179326,-9.40006e-06,1.81362e-09,-1.23836e-13,2076.49,32.2525], Tmin=(868.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = 'C=CO[CH]O(789)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {4,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {6,S}
4  C u1 p0 c0 {1,S} {2,S} {7,S}
5  C u0 p0 c0 {3,D} {8,S} {9,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-150.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.925173,'amu*angstrom^2'), symmetry=1, barrier=(21.2715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925924,'amu*angstrom^2'), symmetry=1, barrier=(21.2888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92581,'amu*angstrom^2'), symmetry=1, barrier=(21.2862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33518,0.0449937,-1.22564e-05,-3.46895e-08,2.21246e-11,-18047,18.1203], Tmin=(100,'K'), Tmax=(901.986,'K')), NASAPolynomial(coeffs=[20.3292,-0.00198739,3.92435e-06,-8.62018e-10,5.8214e-14,-22988.8,-79.9557], Tmin=(901.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJO)"""),
)

species(
    label = '[CH2]COC[O](790)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (11.2052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,541.879,541.887,541.922],'cm^-1')),
        HinderedRotor(inertia=(0.128259,'amu*angstrom^2'), symmetry=1, barrier=(2.94893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141482,'amu*angstrom^2'), symmetry=1, barrier=(2.94888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150063,'amu*angstrom^2'), symmetry=1, barrier=(31.2743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70059,0.0437292,-3.21688e-05,1.17765e-08,-1.72382e-12,1436.16,20.4439], Tmin=(100,'K'), Tmax=(1611.77,'K')), NASAPolynomial(coeffs=[12.8067,0.0161664,-6.51714e-06,1.16627e-09,-7.80489e-14,-2143.9,-38.4377], Tmin=(1611.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.2052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(CJCO)"""),
)

species(
    label = 'C[CH]O[CH]O(791)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {11,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u1 p0 c0 {1,S} {3,S} {9,S}
5  C u1 p0 c0 {1,S} {2,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-58.3656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16954,0.0624009,-8.0762e-05,5.3888e-08,-1.4043e-11,-6917.97,20.4624], Tmin=(100,'K'), Tmax=(946.418,'K')), NASAPolynomial(coeffs=[12.3848,0.0149995,-5.63344e-06,9.6602e-10,-6.32764e-14,-9040.81,-33.0269], Tmin=(946.418,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.3656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + radical(CCsJOCs) + radical(OCJO)"""),
)

species(
    label = '[CH2][CH]OCO(792)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
4  C u1 p0 c0 {1,S} {5,S} {8,S}
5  C u1 p0 c0 {4,S} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-35.7367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.246672,0.0664228,-7.19197e-05,3.54156e-08,-6.26887e-12,-4121.87,26.4103], Tmin=(100,'K'), Tmax=(1689.55,'K')), NASAPolynomial(coeffs=[19.7358,0.000435645,3.24752e-06,-7.87473e-10,5.62528e-14,-8208.17,-72.5843], Tmin=(1689.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.7367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = 'CCO[CH][O](793)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {2,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-11.4237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,180,957.066,957.698],'cm^-1')),
        HinderedRotor(inertia=(0.0975991,'amu*angstrom^2'), symmetry=1, barrier=(2.244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0986766,'amu*angstrom^2'), symmetry=1, barrier=(2.26877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971604,'amu*angstrom^2'), symmetry=1, barrier=(2.23391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24718,0.0287734,-1.15344e-05,1.12071e-09,6.02937e-14,-1418.39,10.2438], Tmin=(100,'K'), Tmax=(2745.94,'K')), NASAPolynomial(coeffs=[38.9552,-0.00824206,1.28738e-06,-1.96021e-10,1.6642e-14,-25585.6,-201.559], Tmin=(2745.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.4237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
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
    E0 = (-40.3057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (381.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (354.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-32.0213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (23.0945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-15.3324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (5.45331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (52.7394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (61.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (106.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (310.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (101.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (118.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (102.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (35.3971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CO[CH]O(784)'],
    products = ['O=CO(77)', 'C2H4(29)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]O(118)', '[CH2]C[O](66)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(17)', '[CH2]O[CH]O(64)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CO[CH]O(784)'],
    products = ['OC1CCO1(785)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CO[CH]O(784)'],
    products = ['C=COCO(788)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CO[CH]O(784)'],
    products = ['CCOC=O(259)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=CO(77)', '[CH2][CH2](81)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.78129e+07,'m^3/(mol*s)'), n=-0.198276, Ea=(89.8628,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3273107057168037, var=7.214547064829838, Tref=1000.0, N=427, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', '[CH2]COC=O(225)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.07,'m^3/(mol*s)'), n=2.43, Ea=(50.0368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][CH]O(212)', 'C2H4(29)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(27.525,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', 'C=CO[CH]O(789)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.68706e+22,'m^3/(mol*s)'), n=-4.56653, Ea=(58.5332,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.23436359975141804, var=0.013647226894766704, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]O(212)', '[CH2][CH2](81)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]COC[O](790)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.41377e+07,'s^-1'), n=1.57383, Ea=(103.529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CO[CH]O(784)'],
    products = ['C[CH]O[CH]O(791)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.7e+13,'s^-1','+|-',2), n=-0.1, Ea=(158.364,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""From training reaction 347 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CO[CH]O(784)'],
    products = ['[CH2][CH]OCO(792)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeO;XH_out] for rate rule [R3H_SS_O;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCO[CH][O](793)'],
    products = ['[CH2]CO[CH]O(784)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #413',
    isomers = [
        '[CH2]CO[CH]O(784)',
    ],
    reactants = [
        ('O=CO(77)', 'C2H4(29)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #413',
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

