species(
    label = '[OH+]=[C-]C1(F)C=C1(8875)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p1 c+1 {6,D} {9,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {3,S} {4,D} {8,S}
6 C u0 p1 c-1 {2,D} {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (354.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,191.399,191.403,191.411,191.424,1278.54,1278.54,1278.54,1278.54,1278.54,3043],'cm^-1')),
        HinderedRotor(inertia=(0.279432,'amu*angstrom^2'), symmetry=1, barrier=(7.26645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91532,0.0531495,-8.99297e-05,8.72817e-08,-3.29663e-11,42725.8,17.8719], Tmin=(100,'K'), Tmax=(817.022,'K')), NASAPolynomial(coeffs=[4.20252,0.0268772,-1.40195e-05,2.75844e-09,-1.92995e-13,42855.2,10.3787], Tmin=(817.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + ring(Cd-Cs(F)-Cd)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC1=CC1(5834)',
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
    label = '[OH+]=[C-]C1=CC1F(11776)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p1 c+1 {6,D} {9,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {3,S} {4,D} {8,S}
6 C u0 p1 c-1 {2,D} {4,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (377.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0075,0.050811,-8.42215e-05,8.24883e-08,-3.16737e-11,45505.7,17.7155], Tmin=(100,'K'), Tmax=(807.502,'K')), NASAPolynomial(coeffs=[3.79325,0.0276047,-1.44382e-05,2.85257e-09,-2.00392e-13,45685.4,12.3811], Tmin=(807.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + ring(Cd-Cs(F)-Cd)"""),
)

species(
    label = '[OH+]=[C-]C1C=C1F(8597)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p1 c+1 {6,D} {9,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {3,S} {4,D} {8,S}
6 C u0 p1 c-1 {2,D} {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (419.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,323,467,575,827,1418,345.509,835.833,835.833,835.833,835.833,835.833,1713.09,1713.09,1713.09,1713.09,1713.09],'cm^-1')),
        HinderedRotor(inertia=(0.124594,'amu*angstrom^2'), symmetry=1, barrier=(3.43206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0813,0.0500455,-8.54988e-05,8.59357e-08,-3.32857e-11,50512.9,17.3765], Tmin=(100,'K'), Tmax=(822.362,'K')), NASAPolynomial(coeffs=[2.85222,0.0288355,-1.49637e-05,2.9368e-09,-2.05115e-13,50976.5,17.3976], Tmin=(822.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = 'F[OH+][C-]=C1C=C1(11777)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p1 c+1 {1,S} {6,S} {9,S}
3 C u0 p0 c0 {4,S} {5,S} {6,D}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {3,S} {4,D} {8,S}
6 C u0 p1 c-1 {2,S} {3,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (681.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,3150,900,1100,612.513,1139.76,1139.76,1139.76,1139.76,1139.76,1139.76,1139.76,1139.76,1139.76,1139.76,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0532497,'amu*angstrom^2'), symmetry=1, barrier=(5.38652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94787,0.0515163,-9.09957e-05,8.45312e-08,-3.00997e-11,81980.4,15.8946], Tmin=(100,'K'), Tmax=(840.333,'K')), NASAPolynomial(coeffs=[6.16814,0.0186648,-9.57385e-06,1.86269e-09,-1.28834e-13,81721.7,-1.05041], Tmin=(840.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CsJ2_singlet-CsH) + ring(methylenecyclopropene)"""),
)

species(
    label = 'FC1=[C-][OH+]C=C1(11778)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p1 c+1 {3,S} {6,S} {9,S}
3 C u0 p0 c0 {2,S} {4,D} {7,S}
4 C u0 p0 c0 {3,D} {5,S} {8,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u0 p1 c-1 {2,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (228.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81413,0.0195258,9.50939e-06,-2.22747e-08,8.50573e-12,27487.7,13.507], Tmin=(100,'K'), Tmax=(1077,'K')), NASAPolynomial(coeffs=[7.98335,0.0158797,-7.07304e-06,1.39785e-09,-1.01624e-13,25472.3,-16.0024], Tmin=(1077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(CsJ2_singlet-CsH) + ring(Cyclopentane)"""),
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
    label = '[OH+]=[C-]C1=C=C1(8879)',
    structure = adjacencyList("""1 O u0 p1 c+1 {5,D} {7,S}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u0 p0 c0 {2,S} {4,D} {6,S}
4 C u0 p0 c0 {2,D} {3,D}
5 C u0 p1 c-1 {1,D} {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (812.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0579,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18542,0.0468254,-8.62326e-05,8.26112e-08,-2.98838e-11,97751.8,12.9449], Tmin=(100,'K'), Tmax=(854.779,'K')), NASAPolynomial(coeffs=[4.96728,0.0184775,-9.58476e-06,1.85e-09,-1.2682e-13,97836.3,3.23649], Tmin=(854.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(812.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + group(CsJ2_singlet-CsH) + ring(Cyclopropadiene)"""),
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
    E0 = (44.3449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (344.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (450.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (520.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (126.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (384.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[OH+]=[C-]C1(F)C=C1(8875)'],
    products = ['CO(13)', 'FC1=CC1(5834)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(702.966,'s^-1'), n=2.72887, Ea=(20.5239,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[OH+]=[C-]C1=CC1F(11776)'],
    products = ['[OH+]=[C-]C1(F)C=C1(8875)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(297.177,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[OH+]=[C-]C1C=C1F(8597)'],
    products = ['[OH+]=[C-]C1(F)C=C1(8875)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(361.847,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[OH+][C-]=C1C=C1(11777)'],
    products = ['[OH+]=[C-]C1(F)C=C1(8875)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(170.518,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[OH+]=[C-]C1(F)C=C1(8875)'],
    products = ['FC1=[C-][OH+]C=C1(11778)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.87873e+12,'s^-1'), n=0.324012, Ea=(103.056,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HF(38)', '[OH+]=[C-]C1=C=C1(8879)'],
    products = ['[OH+]=[C-]C1(F)C=C1(8875)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(184.551,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #3914',
    isomers = [
        '[OH+]=[C-]C1(F)C=C1(8875)',
    ],
    reactants = [
        ('CO(13)', 'FC1=CC1(5834)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3914',
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

