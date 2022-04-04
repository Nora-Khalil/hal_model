species(
    label = 'C=CC[C-]=[OH+](328)',
    structure = adjacencyList("""1  O u0 p1 c+1 {5,D} {11,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  C u0 p1 c-1 {1,D} {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (332.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180,1842.76,1843.16,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(1.0647,'amu*angstrom^2'), symmetry=1, barrier=(24.4794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.370277,'amu*angstrom^2'), symmetry=1, barrier=(24.4811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79945,0.0299711,-1.37988e-05,2.37762e-09,-1.32416e-13,39936.1,12.2409], Tmin=(100,'K'), Tmax=(2977.1,'K')), NASAPolynomial(coeffs=[39.1527,-0.00985374,2.39956e-06,-3.83703e-10,2.6742e-14,15484.7,-202.599], Tmin=(2977.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C=[C-][OH+]C=C(364)',
    structure = adjacencyList("""1  O u0 p1 c+1 {2,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {6,S}
3  C u0 p0 c0 {2,D} {8,S} {9,S}
4  C u0 p0 c0 {5,D} {10,S} {11,S}
5  C u0 p1 c-1 {1,S} {4,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (386.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.346088,'amu*angstrom^2'), symmetry=1, barrier=(7.95724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344831,'amu*angstrom^2'), symmetry=1, barrier=(7.92833,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09266,0.0500565,-8.33543e-05,8.41299e-08,-3.25832e-11,46532.1,16.4693], Tmin=(100,'K'), Tmax=(833.06,'K')), NASAPolynomial(coeffs=[2.23432,0.0308278,-1.53331e-05,2.9676e-09,-2.05719e-13,47152.2,19.6747], Tmin=(833.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CC=[C-][OH2+](365)',
    structure = adjacencyList("""1  O u0 p1 c+1 {5,S} {10,S} {11,S}
2  C u0 p0 c0 {3,S} {4,D} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p1 c-1 {1,S} {3,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (381.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.250563,'amu*angstrom^2'), symmetry=1, barrier=(5.76094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246739,'amu*angstrom^2'), symmetry=1, barrier=(5.67301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20253,0.0442801,-6.10306e-05,5.63446e-08,-2.14928e-11,45989.6,17.2081], Tmin=(100,'K'), Tmax=(786.929,'K')), NASAPolynomial(coeffs=[3.88435,0.0274829,-1.32901e-05,2.58022e-09,-1.80763e-13,45980.3,11.12], Tmin=(786.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C=[C-][OH2+](366)',
    structure = adjacencyList("""1 O u0 p1 c+1 {3,S} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p1 c-1 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (322.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34917,0.0191954,-3.34749e-05,3.95295e-08,-1.69785e-11,38815.2,9.15236], Tmin=(100,'K'), Tmax=(831.117,'K')), NASAPolynomial(coeffs=[0.881183,0.0184289,-9.27096e-06,1.80956e-09,-1.26144e-13,39662.1,23.2297], Tmin=(831.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (73.6726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (252.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (298.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (357.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=CC[C-]=[OH+](328)'],
    products = ['CO(14)', 'C=CC(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(10.6171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C-][OH+]C=C(364)'],
    products = ['C=CC[C-]=[OH+](328)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(135.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC=[C-][OH2+](365)'],
    products = ['C=CC[C-]=[OH+](328)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(150.201,'s^-1'), n=3.11041, Ea=(185.423,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00848789959541425, var=6.520606768146176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R!H->C',), comment="""Estimated from node Root_3R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC[C-]=[OH+](328)'],
    products = ['C2H2(22)', 'C=[C-][OH2+](366)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.8901e+09,'s^-1'), n=0.601596, Ea=(294.517,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.07591793873856297, var=4.052420941329839, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_5R!H->O',), comment="""Estimated from node Root_1R!H->C_2R!H->C_5R!H->O
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #242',
    isomers = [
        'C=CC[C-]=[OH+](328)',
    ],
    reactants = [
        ('CO(14)', 'C=CC(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #242',
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

