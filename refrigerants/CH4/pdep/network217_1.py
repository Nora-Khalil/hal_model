species(
    label = 'O=C1OOC1=O(267)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
"""),
    E0 = (-463.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551.294,551.294,551.294,551.294,551.294,551.294,551.295,551.295,551.295,551.295,551.295,551.295],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03638,-0.00926841,0.00013209,-1.7716e-07,6.94151e-11,-55702.5,13.9659], Tmin=(100,'K'), Tmax=(948.188,'K')), NASAPolynomial(coeffs=[20.7701,-0.00621505,4.08134e-06,-5.51469e-10,1.57863e-14,-62565.7,-89.1028], Tmin=(948.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-463.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)O2s) + ring(Cyclobutane)"""),
)

species(
    label = 'CO2(15)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=C=C1OOO1(1248)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (157.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,564.009,564.009,564.009,564.009,564.009,564.009,564.009,564.009,564.009],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67981,-0.000742596,0.000108526,-1.64065e-07,6.91735e-11,19065.6,17.4177], Tmin=(100,'K'), Tmax=(918.039,'K')), NASAPolynomial(coeffs=[24.5557,-0.0188155,1.18471e-05,-2.20729e-09,1.38071e-13,11794,-103.978], Tmin=(918.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + missing(Cdd-CdO2d) + ring(Cyclobutane)"""),
)

species(
    label = 'O1O[C]2OO[C]12(1249)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u1 p0 c0 {1,S} {3,S} {6,S}
6 C u1 p0 c0 {2,S} {4,S} {5,S}
"""),
    E0 = (413.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([646.642,646.778,646.808,646.937,646.945,646.949,646.983,646.983,646.986,647.002,647.009,647.096],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.78356,0.0181054,1.40303e-05,-2.84117e-08,1.02643e-11,49779.5,14.2872], Tmin=(100,'K'), Tmax=(1142.54,'K')), NASAPolynomial(coeffs=[10.8915,0.0123316,-8.07558e-06,1.80873e-09,-1.38479e-13,46450.9,-32.3679], Tmin=(1142.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + polycyclic(s2_4_4_ane) + radical(Cs_P) + radical(Cs_P)"""),
)

species(
    label = 'O=[C]OO[C]=O(266)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {3,D}
6 C u1 p0 c0 {2,S} {4,D}
"""),
    E0 = (-214.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1850,1860,440,470,900,1000],'cm^-1')),
        HinderedRotor(inertia=(1.89499,'amu*angstrom^2'), symmetry=1, barrier=(43.5696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88973,'amu*angstrom^2'), symmetry=1, barrier=(43.4486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3661.19,'J/mol'), sigma=(5.65648,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.87 K, Pc=45.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06027,0.0235154,3.48158e-05,-7.80732e-08,3.54793e-11,-25689.8,18.2992], Tmin=(100,'K'), Tmax=(937.227,'K')), NASAPolynomial(coeffs=[20.392,-0.00598464,4.02568e-06,-6.8606e-10,3.62432e-14,-31266.6,-80.3722], Tmin=(937.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(120.56,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(Cds-OdOsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical((O)CJOC)"""),
)

species(
    label = '[O]OC(=O)[C]=O(1250)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 O u1 p2 c0 {1,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,D} {6,S}
6 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-52.7583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1855,455,950,180,180,180,3584.97],'cm^-1')),
        HinderedRotor(inertia=(0.585677,'amu*angstrom^2'), symmetry=1, barrier=(13.4659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585366,'amu*angstrom^2'), symmetry=1, barrier=(13.4587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75401,0.0533816,-9.73729e-05,8.01671e-08,-2.43577e-11,-6268.14,18.649], Tmin=(100,'K'), Tmax=(954.724,'K')), NASAPolynomial(coeffs=[10.674,0.00416036,-1.42317e-06,1.67584e-10,-5.38176e-15,-7431.35,-21.1436], Tmin=(954.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.7583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C(=O)OOJ) + radical(OC=OCJ=O)"""),
)

species(
    label = '[O]C(=O)C([O])=O(264)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {5,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
"""),
    E0 = (-217.699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4825.3,'J/mol'), sigma=(6.37974,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=753.70 K, Pc=42.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57497,0.0390004,-4.27579e-05,-1.76496e-08,4.81105e-11,-26139.8,17.0519], Tmin=(100,'K'), Tmax=(457.299,'K')), NASAPolynomial(coeffs=[7.03639,0.0143462,-9.02393e-06,1.88718e-09,-1.36015e-13,-26698.1,-2.62398], Tmin=(457.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ) + radical(C=OC=OOJ)"""),
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
    E0 = (-336.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (428.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (540.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-79.1704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (82.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-82.5444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CO2(15)', 'CO2(15)'],
    products = ['O=C1OOC1=O(267)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.46542e-32,'m^3/(mol*s)'), n=10.9581, Ea=(342.618,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1COCSCdCdd->Cd',), comment="""Estimated from node Root_N-1COCSCdCdd->Cd
Multiplied by reaction path degeneracy 2.0
Ea raised from 335.0 to 342.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C=C1OOO1(1248)'],
    products = ['O=C1OOC1=O(267)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.19768e+62,'s^-1'), n=-13.6281, Ea=(144.126,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O1O[C]2OO[C]12(1249)'],
    products = ['O=C1OOC1=O(267)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]OO[C]=O(266)'],
    products = ['O=C1OOC1=O(267)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(=O)[C]=O(1250)'],
    products = ['O=C1OOC1=O(267)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(=O)C([O])=O(264)'],
    products = ['O=C1OOC1=O(267)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.48e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 4.0
family: Birad_recombination"""),
)

network(
    label = 'PDepNetwork #217',
    isomers = [
        'O=C1OOC1=O(267)',
    ],
    reactants = [
        ('CO2(15)', 'CO2(15)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #217',
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

