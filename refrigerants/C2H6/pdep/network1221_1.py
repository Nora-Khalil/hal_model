species(
    label = 'CC#CC[C-]=[OH+](2419)',
    structure = adjacencyList("""1  O u0 p1 c+1 {6,D} {12,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,T}
5  C u0 p0 c0 {3,S} {4,T}
6  C u0 p1 c-1 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (455.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65672,0.0590497,-8.80964e-05,8.06735e-08,-2.85655e-11,54919.6,21.1378], Tmin=(100,'K'), Tmax=(887.692,'K')), NASAPolynomial(coeffs=[3.38791,0.0333367,-1.43797e-05,2.58026e-09,-1.70393e-13,55318,16.9671], Tmin=(887.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CO(14)',
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
    label = 'C=C=CC(1692)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.759585,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74634,0.0218191,8.22302e-06,-2.14762e-08,8.55596e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.61,'K')), NASAPolynomial(coeffs=[6.82084,0.0192337,-7.45616e-06,1.36535e-09,-9.53184e-14,16028,-10.4336], Tmin=(1025.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C#CC[C-]=[OH+](1426)',
    structure = adjacencyList("""1 O u0 p1 c+1 {3,D} {8,S}
2 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3 C u0 p1 c-1 {1,D} {2,S}
4 C u0 p0 c0 {2,S} {5,T}
5 C u0 p0 c0 {4,T} {9,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (498.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100,180,180,180,1489.68,1492,1497.37],'cm^-1')),
        HinderedRotor(inertia=(0.245775,'amu*angstrom^2'), symmetry=1, barrier=(5.65085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244552,'amu*angstrom^2'), symmetry=1, barrier=(5.62272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0738,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2836.07,'J/mol'), sigma=(5.784e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0506,0.0511631,-8.98119e-05,8.66312e-08,-3.12589e-11,59977.9,16.3028], Tmin=(100,'K'), Tmax=(883.874,'K')), NASAPolynomial(coeffs=[3.42171,0.0250228,-1.16182e-05,2.13549e-09,-1.4207e-13,60514.2,14.2623], Tmin=(883.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(CsJ2_singlet-CsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C=[C-][OH+]C#CC(2706)',
    structure = adjacencyList("""1  O u0 p1 c+1 {5,S} {6,S} {10,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {5,D} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,T}
5  C u0 p1 c-1 {1,S} {3,D}
6  C u0 p0 c0 {1,S} {4,T}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (502.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2100,2250,500,550,180,180,180,180,1170.57],'cm^-1')),
        HinderedRotor(inertia=(0.216454,'amu*angstrom^2'), symmetry=1, barrier=(4.97671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217972,'amu*angstrom^2'), symmetry=1, barrier=(5.01161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219173,'amu*angstrom^2'), symmetry=1, barrier=(5.03923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70457,0.0584504,-8.94779e-05,8.29193e-08,-2.9405e-11,60571.1,29.5474], Tmin=(100,'K'), Tmax=(894.492,'K')), NASAPolynomial(coeffs=[3.09583,0.0328571,-1.40748e-05,2.51009e-09,-1.64794e-13,61097.2,27.3224], Tmin=(894.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CtHHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(CsJ2_singlet-CsH) + group(Ct-CtCs)"""),
)

species(
    label = 'CC#CC=[C-][OH2+](2707)',
    structure = adjacencyList("""1  O u0 p1 c+1 {5,S} {11,S} {12,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {5,D} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {6,T}
5  C u0 p1 c-1 {1,S} {3,D}
6  C u0 p0 c0 {3,S} {4,T}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (516.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2100,2250,500,550,180,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.112204,'amu*angstrom^2'), symmetry=1, barrier=(2.5798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113206,'amu*angstrom^2'), symmetry=1, barrier=(2.60283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888312,'amu*angstrom^2'), symmetry=1, barrier=(20.424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93584,0.0519948,-7.55972e-05,7.19528e-08,-2.73297e-11,62186.6,20.4801], Tmin=(100,'K'), Tmax=(824.967,'K')), NASAPolynomial(coeffs=[2.96537,0.0330227,-1.56815e-05,2.99221e-09,-2.06623e-13,62492.5,18.5946], Tmin=(824.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(CsJ2_singlet-CsH) + group(Ct-Ct(Cds-Cds))"""),
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
    E0 = (67.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (516.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (240.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (296.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC#CC[C-]=[OH+](2419)'],
    products = ['CO(14)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(12.5308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', 'C#CC[C-]=[OH+](1426)'],
    products = ['CC#CC[C-]=[OH+](2419)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.57523e+09,'m^3/(mol*s)'), n=-0.891994, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_Ext-7R!H-R_N-4CbCdCt->Cd',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-6R!H-R_N-Sp-6R!H-4CbCdCt_N-Sp-7R!H=6R!H_Ext-7R!H-R_N-4CbCdCt->Cd"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C-][OH+]C#CC(2706)'],
    products = ['CC#CC[C-]=[OH+](2419)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(138.508,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CC#CC=[C-][OH2+](2707)'],
    products = ['CC#CC[C-]=[OH+](2419)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(150.201,'s^-1'), n=3.11041, Ea=(180.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00848789959541425, var=6.520606768146176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R!H->C',), comment="""Estimated from node Root_3R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1221',
    isomers = [
        'CC#CC[C-]=[OH+](2419)',
    ],
    reactants = [
        ('CO(14)', 'C=C=CC(1692)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1221',
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

