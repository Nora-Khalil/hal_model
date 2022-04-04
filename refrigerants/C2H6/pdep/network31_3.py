species(
    label = '[CH]=CC[CH2](65)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u1 p0 c0 {1,S} {8,S} {9,S}
4  C u1 p0 c0 {2,D} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (435.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.207738,'amu*angstrom^2'), symmetry=1, barrier=(4.77631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204754,'amu*angstrom^2'), symmetry=1, barrier=(4.7077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52057,0.0279428,-9.18493e-06,-4.19504e-09,2.67709e-12,52484.1,16.9976], Tmin=(100,'K'), Tmax=(1112.3,'K')), NASAPolynomial(coeffs=[7.71989,0.0177227,-6.83491e-06,1.24846e-09,-8.64342e-14,50803,-10.9969], Tmin=(1112.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=C(61)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,D} {5,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {7,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80601,0.0102582,6.17268e-05,-9.01653e-08,3.59121e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.045,'K')), NASAPolynomial(coeffs=[12.4693,0.0100555,-2.41212e-06,4.5709e-10,-3.93171e-14,8010.8,-43.6372], Tmin=(946.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C1=CCC1(69)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (144.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,204.182,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,3484.69],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3170.39,'J/mol'), sigma=(5.50295,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.21 K, Pc=43.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39615,-0.00343641,9.09434e-05,-1.1327e-07,4.24092e-11,17412.4,10.7746], Tmin=(100,'K'), Tmax=(950.248,'K')), NASAPolynomial(coeffs=[9.42764,0.013964,-4.06869e-06,7.74979e-10,-6.20376e-14,14334.2,-28.1812], Tmin=(950.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C2H2(22)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'C2H3(28)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=[CH](58)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,D} {3,S}
2 C u1 p0 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1227.31,2941.71],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88058,-0.000843434,2.04094e-05,-2.51944e-08,9.47575e-12,71059.3,4.72115], Tmin=(100,'K'), Tmax=(935.367,'K')), NASAPolynomial(coeffs=[4.92141,0.00358272,-9.24433e-07,1.57271e-10,-1.19539e-14,70476.2,-2.30656], Tmin=(935.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=CC=C(66)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,D} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 C u1 p0 c0 {2,D} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0824,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64256,0.0163337,3.86225e-05,-6.71377e-08,2.83603e-11,41729.6,13.282], Tmin=(100,'K'), Tmax=(937.724,'K')), NASAPolynomial(coeffs=[12.9705,0.00669127,-1.00071e-06,1.67602e-10,-1.71437e-14,38279.7,-43.9475], Tmin=(937.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH2](39)',
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
    label = 'C#CC[CH2](70)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {7,S} {8,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {3,T} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (357.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.363298,'amu*angstrom^2'), symmetry=1, barrier=(8.35294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0523376,'amu*angstrom^2'), symmetry=1, barrier=(27.2367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0824,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66701,0.0244369,-3.88217e-06,-1.05367e-08,5.37535e-12,43000.9,14.3099], Tmin=(100,'K'), Tmax=(1010.88,'K')), NASAPolynomial(coeffs=[7.82292,0.0147589,-5.43381e-06,9.80744e-10,-6.83221e-14,41410.6,-13.3303], Tmin=(1010.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ)"""),
)

species(
    label = '[CH]C=CC(71)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  C u2 p0 c0 {3,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (342.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,213.499,214.183,214.544],'cm^-1')),
        HinderedRotor(inertia=(1.57513,'amu*angstrom^2'), symmetry=1, barrier=(51.0829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56845,'amu*angstrom^2'), symmetry=1, barrier=(51.0769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52,0.0276811,9.44242e-07,-1.25014e-08,4.67165e-12,41238.4,14.8113], Tmin=(100,'K'), Tmax=(1157.78,'K')), NASAPolynomial(coeffs=[6.04216,0.0260074,-1.04844e-05,1.90883e-09,-1.30602e-13,39719.5,-5.73484], Tmin=(1157.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C[C]=C(72)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u1 p0 c0 {1,S} {7,S} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (426.648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,2857.61],'cm^-1')),
        HinderedRotor(inertia=(0.288735,'amu*angstrom^2'), symmetry=1, barrier=(6.63859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842868,'amu*angstrom^2'), symmetry=1, barrier=(19.3792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5849,0.0285817,-1.50813e-05,3.88347e-09,-3.99393e-13,51367,16.8487], Tmin=(100,'K'), Tmax=(2186,'K')), NASAPolynomial(coeffs=[10.4562,0.0141787,-5.19829e-06,8.69449e-10,-5.47005e-14,47925.7,-27.2818], Tmin=(2186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC(73)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {4,D}
4  C u1 p0 c0 {3,D} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (468.498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.366722,'amu*angstrom^2'), symmetry=1, barrier=(8.43167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366981,'amu*angstrom^2'), symmetry=1, barrier=(8.43761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59813,0.0301253,-1.77975e-05,5.47286e-09,-7.04425e-13,56398,15.4011], Tmin=(100,'K'), Tmax=(1698.28,'K')), NASAPolynomial(coeffs=[7.75978,0.0179682,-7.05993e-06,1.25786e-09,-8.39556e-14,54644.8,-12.2346], Tmin=(1698.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C[CH2](74)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,D} {3,S} {5,S}
2  C u0 p0 c0 {1,D} {4,S} {6,S}
3  C u1 p0 c0 {1,S} {7,S} {10,S}
4  C u1 p0 c0 {2,S} {8,S} {9,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.209873,'amu*angstrom^2'), symmetry=1, barrier=(25.2262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778586,'amu*angstrom^2'), symmetry=1, barrier=(93.3137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56317,0.0223429,1.87066e-05,-3.93097e-08,1.63982e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.265,'K')), NASAPolynomial(coeffs=[9.82996,0.0151966,-5.22271e-06,9.67654e-10,-7.07858e-14,30607.7,-26.985], Tmin=(974.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C=C(1095)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {4,D} {8,S} {9,S}
4 C u1 p0 c0 {1,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (303.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.39164,'amu*angstrom^2'), symmetry=1, barrier=(31.9966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0824,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67736,0.0217333,1.0387e-05,-2.87121e-08,1.2538e-11,36547.1,13.0516], Tmin=(100,'K'), Tmax=(966.654,'K')), NASAPolynomial(coeffs=[8.97812,0.0135587,-4.70087e-06,8.4736e-10,-6.05066e-14,34492.7,-21.4574], Tmin=(966.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C[C]C(1830)',
    structure = adjacencyList("""1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,D} {4,S} {8,S}
3  C u0 p0 c0 {2,D} {9,S} {10,S}
4  C u0 p1 c0 {1,S} {2,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (411.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0841,'amu*angstrom^2'), symmetry=1, barrier=(47.9177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08658,'amu*angstrom^2'), symmetry=1, barrier=(47.9745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68372,0.0299287,-1.41069e-05,2.85821e-09,-2.07263e-13,49479.1,23.4398], Tmin=(100,'K'), Tmax=(2271.38,'K')), NASAPolynomial(coeffs=[13.3021,0.0141313,-5.59086e-06,9.21172e-10,-5.59724e-14,43906.9,-38.1466], Tmin=(2271.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC=C-2(1831)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,D} {7,S}
3  C u0 p0 c0 {2,D} {8,S} {9,S}
4  C u0 p1 c0 {1,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (433.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,338.785,338.803,2090.38,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0834505,'amu*angstrom^2'), symmetry=1, barrier=(6.79651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329161,'amu*angstrom^2'), symmetry=1, barrier=(26.8194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67988,0.0290194,-1.49085e-05,3.4546e-09,-3.0962e-13,52178.5,14.2992], Tmin=(100,'K'), Tmax=(2560.9,'K')), NASAPolynomial(coeffs=[14.0172,0.0113115,-4.53664e-06,7.54625e-10,-4.60496e-14,46371.6,-51.0584], Tmin=(2560.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C1=CC1(798)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 C u0 p0 c0 {1,S} {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (271.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,918.87,918.871,918.872,918.872,918.872,918.873,3181.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.6352,-0.00323151,6.10389e-05,-7.68404e-08,2.92237e-11,32671.6,7.3571], Tmin=(100,'K'), Tmax=(938.64,'K')), NASAPolynomial(coeffs=[7.97924,0.00707126,-1.47314e-06,2.63802e-10,-2.31982e-14,30586.7,-20.087], Tmin=(938.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""cC3H4""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]1C[CH]C1(2182)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
3  C u1 p0 c0 {1,S} {2,S} {9,S}
4  C u1 p0 c0 {1,S} {2,S} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (389.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,180,388.393,1395.57,1395.58,1395.59,1395.6,1395.6,1395.6,1395.6,1395.61,1395.61,1802.73],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61382,-0.00133354,6.67738e-05,-7.44503e-08,2.51141e-11,46874.5,14.5614], Tmin=(100,'K'), Tmax=(1012.31,'K')), NASAPolynomial(coeffs=[3.96904,0.0238095,-9.81749e-06,1.89418e-09,-1.37336e-13,45442.4,6.12497], Tmin=(1012.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1[CH]CC1(2183)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {4,S} {9,S}
4  C u1 p0 c0 {2,S} {3,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (389.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,180,388.393,1395.57,1395.58,1395.59,1395.6,1395.6,1395.6,1395.6,1395.61,1395.61,1802.73],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61382,-0.00133354,6.67738e-05,-7.44503e-08,2.51141e-11,46874.5,15.2545], Tmin=(100,'K'), Tmax=(1012.31,'K')), NASAPolynomial(coeffs=[3.96904,0.0238095,-9.81749e-06,1.89418e-09,-1.37336e-13,45442.4,6.81812], Tmin=(1012.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1C=CC1(1092)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {4,S} {8,S}
3 C u0 p0 c0 {1,S} {4,D} {7,S}
4 C u0 p0 c0 {2,S} {3,D} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (309.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,435.348,1132.79,1133.03,1133.39,1133.5,1133.72,1134.19,1134.34,1134.39,1134.44,1135.48],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0824,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59806,-0.00870546,0.000100053,-1.21959e-07,4.54681e-11,37195.9,8.94337], Tmin=(100,'K'), Tmax=(948.075,'K')), NASAPolynomial(coeffs=[9.48648,0.0114517,-3.03707e-06,5.96472e-10,-5.06745e-14,34056.9,-29.8169], Tmin=(948.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[C]1=CCC1(2184)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3 C u0 p0 c0 {1,S} {4,D} {9,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (396.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,180,589.984,1114.93,1114.93,1114.94,1114.95,1114.95,1114.95,1114.95,1114.96,3358.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0824,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36702,0.00111968,6.62926e-05,-8.41431e-08,3.14826e-11,47776.9,12.046], Tmin=(100,'K'), Tmax=(956.527,'K')), NASAPolynomial(coeffs=[8.22875,0.0134799,-4.35533e-06,8.26088e-10,-6.35385e-14,45351.3,-19.0104], Tmin=(956.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[C]1CCC1(2185)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p1 c0 {2,S} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (442.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,875.825,875.825,875.825,875.825,875.825,875.825,875.825,875.826,875.826,875.826,875.826,875.826],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16729,0.00609402,6.26578e-05,-7.90276e-08,2.87147e-11,53265.8,21.4593], Tmin=(100,'K'), Tmax=(980.338,'K')), NASAPolynomial(coeffs=[7.14685,0.0210755,-8.03279e-06,1.52847e-09,-1.12012e-13,50985.3,-5.31204], Tmin=(980.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsJ2_singlet-CsH) + ring(Cyclobutane)"""),
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
    E0 = (109.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (117.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (172.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (315.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (231.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (234.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (291.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (582.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (216.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (215.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (294.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (301.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-107.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (123.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (41.0435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (246.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (193.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (231.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (110.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (119.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (482.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-43.3315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (86.3661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (126.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (194.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (282.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (156.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC[CH2](65)'],
    products = ['C2H2(22)', 'C2H4(29)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CC[CH2](65)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CC[CH2](65)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[CH](58)', 'C2H4(29)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(309.517,'m^3/(mol*s)'), n=1.24899, Ea=(8.66466,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.02821839893876885, var=2.360387306148194, Tref=1000.0, N=229, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(6)', '[CH]=CC=C(66)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000216213,'m^3/(mol*s)'), n=3.29826, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.02158166639280058, var=0.8466995369935014, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_N-Sp-5R!H-4CClOS_Sp-4CClOS-1CCClOOSS_Sp-2CS=1CCOSS_Sp-5R!H=4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_N-Sp-5R!H-4CClOS_Sp-4CClOS-1CCClOOSS_Sp-2CS=1CCOSS_Sp-5R!H=4CClOS"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H2(22)', '[CH2][CH2](39)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0272,'m^3/(mol*s)'), n=2.41, Ea=(24.5255,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_N-1CN->N_N-Sp-2CNO=1C_2CNO->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_N-1CN->N_N-Sp-2CNO=1C_2CNO->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', 'C#CC[CH2](70)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.53376e+24,'m^3/(mol*s)'), n=-5.13718, Ea=(48.7799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_5R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_5R!H->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[CH](58)', '[CH2][CH2](39)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.83229e+09,'m^3/(mol*s)'), n=-0.619182, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.003931094435416637, var=3.092035168239547, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_N-3BrCO->Br_N-3CO->O_Sp-3C-1C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_N-3BrCO->Br_N-3CO->O_Sp-3C-1C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]C=CC(71)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.18e+09,'s^-1'), n=1.22, Ea=(199.995,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 157 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC[CH2](65)'],
    products = ['[CH2]C[C]=C(72)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC(73)'],
    products = ['[CH]=CC[CH2](65)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC[CH2](65)'],
    products = ['[CH2]C=C[CH2](74)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=CC=C(61)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH2_1;CdH2_2]
Euclidian distance = 1.4142135623730951
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C[C]=C(72)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C=CC(71)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction1',
    reactants = ['C2H3(28)', 'C2H3(28)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.615e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Matched reaction 89 C2H3 + C2H3 <=> C4H6-4 in R_Recombination/training
This reaction matched rate rule [Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_N-3BrCO->Br_N-3CO->O_N-Sp-3C-1C_N-Sp-3C#1C_Ext-2CF-R]
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(6)', 'C=[C]C=C(1095)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(1.003,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 49 C4H5 + H <=> C4H6-2 in R_Recombination/training
This reaction matched rate rule [Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_N-Sp-4R!H-2C]
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(6)', '[CH]=CC=C(66)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.83712e+08,'m^3/(mol*s)'), n=-0.251115, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.04900385874705703, var=3.2086570313419123, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C
Ea raised from -5.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[C]C(1830)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.99999e+12,'s^-1'), n=8.2394e-08, Ea=(25.5617,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]CC=C-2(1831)'],
    products = ['C=CC=C(61)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.40886e+17,'s^-1'), n=-1.43042, Ea=(11.9408,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH2(S)(24)', 'C1=CC1(798)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(117.966,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=C[CH2](74)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]1C[CH]C1(2182)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.0412e+11,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 8.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]1[CH]CC1(2183)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(6)', '[CH]1C=CC1(1092)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.40772e+06,'m^3/(mol*s)'), n=0.19941, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-1.8120480258864409, var=13.318773824843154, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(6)', '[C]1=CCC1(2184)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.70386e+06,'m^3/(mol*s)'), n=0.19941, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-1.8120480258864409, var=13.318773824843154, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]1CCC1(2185)'],
    products = ['C1=CCC1(69)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.61387e+21,'s^-1'), n=-2.18552, Ea=(39.9294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5822770821690705, var=6.058307656543806, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C
Multiplied by reaction path degeneracy 4.0"""),
)

network(
    label = 'PDepNetwork #31',
    isomers = [
        '[CH]=CC[CH2](65)',
        'C=CC=C(61)',
        'C1=CCC1(69)',
    ],
    reactants = [
        ('C2H2(22)', 'C2H4(29)'),
        ('C2H3(28)', 'C2H3(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #31',
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

