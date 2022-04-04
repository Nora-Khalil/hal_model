species(
    label = 'C#CC(F)[C-]=[OH+](7376)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p1 c+1 {4,D} {8,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p1 c-1 {2,D} {3,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (299.171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100,180,180,180,180,703.072,703.172],'cm^-1')),
        HinderedRotor(inertia=(0.147638,'amu*angstrom^2'), symmetry=1, barrier=(3.3945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145301,'amu*angstrom^2'), symmetry=1, barrier=(3.34075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5989,0.0638016,-0.000121595,1.19288e-07,-4.34384e-11,36057.9,18.2745], Tmin=(100,'K'), Tmax=(871.617,'K')), NASAPolynomial(coeffs=[3.96245,0.0270185,-1.36586e-05,2.5919e-09,-1.75235e-13,36631.1,12.8481], Tmin=(871.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFH) + group(CsJ2_singlet-CsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
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
    label = 'C=C=CF(6101)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (9.45959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1537.31],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2831,'J/mol'), sigma=(4.74838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=442.20 K, Pc=60 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96414,0.0021816,6.68618e-05,-1.27032e-07,7.57672e-11,1138.94,8.06307], Tmin=(10,'K'), Tmax=(518.59,'K')), NASAPolynomial(coeffs=[2.17835,0.0231528,-1.46137e-05,4.46872e-09,-5.27221e-13,1227.38,14.5728], Tmin=(518.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(9.45959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C#CC=[C-][OH+]F(8080)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p1 c+1 {1,S} {4,S} {8,S}
3 C u0 p0 c0 {4,D} {5,S} {7,S}
4 C u0 p1 c-1 {2,S} {3,D}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (558.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.722144,'amu*angstrom^2'), symmetry=1, barrier=(19.0045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271971,'amu*angstrom^2'), symmetry=1, barrier=(6.25314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31386,0.0443151,-7.81411e-05,7.9161e-08,-3.06492e-11,67245.6,15.7007], Tmin=(100,'K'), Tmax=(824.276,'K')), NASAPolynomial(coeffs=[3.02483,0.0246623,-1.28919e-05,2.54059e-09,-1.77719e-13,67678.8,15.747], Tmin=(824.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[OH+][C-]=CF(8081)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p1 c+1 {4,S} {5,S} {7,S}
3 C u0 p0 c0 {1,S} {4,D} {8,S}
4 C u0 p1 c-1 {2,S} {3,D}
5 C u0 p0 c0 {2,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (353.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.47315,'amu*angstrom^2'), symmetry=1, barrier=(56.8626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523267,'amu*angstrom^2'), symmetry=1, barrier=(12.0309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68064,0.0619236,-0.000118577,1.15196e-07,-4.1146e-11,42575.3,27.7595], Tmin=(100,'K'), Tmax=(892.507,'K')), NASAPolynomial(coeffs=[3.96455,0.0250975,-1.19959e-05,2.20357e-09,-1.45293e-13,43226.7,22.9337], Tmin=(892.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CdCFH) + group(CsJ2_singlet-CsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(F)=[C-][OH2+](8082)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p1 c+1 {4,S} {7,S} {8,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p1 c-1 {2,S} {3,D}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (356.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.800168,'amu*angstrom^2'), symmetry=1, barrier=(18.3974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798062,'amu*angstrom^2'), symmetry=1, barrier=(18.349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83429,0.0554855,-0.000100226,9.82331e-08,-3.71461e-11,42892,17.3805], Tmin=(100,'K'), Tmax=(809.365,'K')), NASAPolynomial(coeffs=[4.89804,0.024537,-1.35738e-05,2.72855e-09,-1.92874e-13,42913.8,6.446], Tmin=(809.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CdCCF) + group(CsJ2_singlet-CsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
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
    E0 = (11.6784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (451.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (186.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (235.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC(F)[C-]=[OH+](7376)'],
    products = ['CO(13)', 'C=C=CF(6101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(15.0313,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC=[C-][OH+]F(8080)'],
    products = ['C#CC(F)[C-]=[OH+](7376)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(195.741,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[OH+][C-]=CF(8081)'],
    products = ['C#CC(F)[C-]=[OH+](7376)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(135.307,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC(F)=[C-][OH2+](8082)'],
    products = ['C#CC(F)[C-]=[OH+](7376)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(150.201,'s^-1'), n=3.11041, Ea=(182.371,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00848789959541425, var=6.520606768146176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R!H->C',), comment="""Estimated from node Root_3R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #2880',
    isomers = [
        'C#CC(F)[C-]=[OH+](7376)',
    ],
    reactants = [
        ('CO(13)', 'C=C=CF(6101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2880',
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

