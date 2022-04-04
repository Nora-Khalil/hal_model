species(
    label = 'C=COC(=O)O(1710)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {6,D} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-531.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44327,0.0416364,-3.00461e-06,-3.73997e-08,2.00421e-11,-63829.3,19.2275], Tmin=(100,'K'), Tmax=(956.196,'K')), NASAPolynomial(coeffs=[18.8942,0.00442612,-7.78712e-07,1.94136e-10,-2.16593e-14,-68802.8,-72.7376], Tmin=(956.196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-531.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'CO2(14)',
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
    label = 'CH3CHO(36)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.19,1305.32,1306.47,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136162,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C(=O)O(309)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-383.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,676.115,680.126,680.289,680.694,683.906],'cm^-1')),
        HinderedRotor(inertia=(0.00541205,'amu*angstrom^2'), symmetry=1, barrier=(1.79464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.87,'J/mol'), sigma=(5.42485,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.51 K, Pc=51.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97086,0.00159581,4.23066e-05,-7.29624e-08,3.78807e-11,-46086.4,8.41125], Tmin=(10,'K'), Tmax=(620.938,'K')), NASAPolynomial(coeffs=[2.60714,0.0173834,-1.27477e-05,4.30843e-09,-5.42128e-13,-46052.1,13.2535], Tmin=(620.938,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-383.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(DO)O""", comment="""Thermo library: 2-BTP_G4"""),
)

species(
    label = 'C2H3(29)',
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
    label = 'O=[C]O(145)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (-192.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1244.76,1684.94],'cm^-1')),
        HinderedRotor(inertia=(0.00103939,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-192.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2CHO(35)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,531.349,531.873,974.177,1641.05,1641.53],'cm^-1')),
        HinderedRotor(inertia=(0.00117252,'amu*angstrom^2'), symmetry=1, barrier=(2.241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CCC(=O)O(1703)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {2,D} {4,S}
6  C u0 p0 c0 {3,D} {4,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
"""),
    E0 = (-538.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27941,0.0361301,-2.28817e-05,6.91437e-09,-8.42719e-13,-64708.6,19.8284], Tmin=(100,'K'), Tmax=(1839.49,'K')), NASAPolynomial(coeffs=[10.8445,0.0175053,-7.69451e-06,1.41029e-09,-9.46835e-14,-67859.7,-26.7138], Tmin=(1839.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'O[C]1O[CH]CO1(2640)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {6,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u1 p0 c0 {2,S} {4,S} {9,S}
6  C u1 p0 c0 {1,S} {2,S} {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-130.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,354.75,1255.83,1255.83,1255.83,1255.83,1255.83,1255.83,1255.83,1255.83,1255.83,1255.83,1255.83,3343.9,3343.9],'cm^-1')),
        HinderedRotor(inertia=(0.0042153,'amu*angstrom^2'), symmetry=1, barrier=(0.44808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41912,0.046475,-5.73759e-05,3.50815e-08,-7.5827e-12,-15536,19.4074], Tmin=(100,'K'), Tmax=(1438.94,'K')), NASAPolynomial(coeffs=[10.3575,0.00500352,3.18529e-06,-1.00578e-09,8.20232e-14,-16387.3,-20.9873], Tmin=(1438.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-OsOsOsH) + ring(1,3-Dioxolane) + radical(CCsJOCs) + radical(Cs_P)"""),
)

species(
    label = 'C=COC([O])[O](2641)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u1 p2 c0 {4,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,D} {8,S}
6  C u0 p0 c0 {5,D} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-74.8439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,988.982,989.116,989.527,989.548],'cm^-1')),
        HinderedRotor(inertia=(0.194258,'amu*angstrom^2'), symmetry=1, barrier=(4.46638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194339,'amu*angstrom^2'), symmetry=1, barrier=(4.46824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33609,0.049465,-9.40461e-05,1.13091e-07,-4.9547e-11,-8954.39,18.7782], Tmin=(100,'K'), Tmax=(810.618,'K')), NASAPolynomial(coeffs=[-3.84758,0.0461501,-2.53147e-05,5.08362e-09,-3.59858e-13,-6840.43,54.168], Tmin=(810.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.8439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = 'C=[C]OC([O])O(2642)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
5  C u0 p0 c0 {6,D} {8,S} {9,S}
6  C u1 p0 c0 {1,S} {5,D}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-62.4976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1685,370,180,180,880.7,880.774],'cm^-1')),
        HinderedRotor(inertia=(0.0063816,'amu*angstrom^2'), symmetry=1, barrier=(3.51357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152812,'amu*angstrom^2'), symmetry=1, barrier=(3.51345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152703,'amu*angstrom^2'), symmetry=1, barrier=(3.51093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.46781,0.0275717,-1.29095e-05,1.83813e-09,-4.22407e-14,-7573.65,14.2471], Tmin=(100,'K'), Tmax=(2884.26,'K')), NASAPolynomial(coeffs=[47.1282,-0.018917,4.67618e-06,-7.03058e-10,4.5965e-14,-37454.2,-245.892], Tmin=(2884.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.4976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCOJ) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]O[C](O)O(2643)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {6,D} {7,S} {8,S}
6  C u1 p0 c0 {1,S} {5,D}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-84.6492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,2950,3100,1380,975,1025,1650,1685,370,180,751.488,4000],'cm^-1')),
        HinderedRotor(inertia=(0.182474,'amu*angstrom^2'), symmetry=1, barrier=(4.19543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697052,'amu*angstrom^2'), symmetry=1, barrier=(16.0266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.040005,'amu*angstrom^2'), symmetry=1, barrier=(16.0272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697064,'amu*angstrom^2'), symmetry=1, barrier=(16.0269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61677,0.055704,-7.42801e-05,5.23576e-08,-1.47816e-11,-10098,24.2456], Tmin=(100,'K'), Tmax=(864.165,'K')), NASAPolynomial(coeffs=[9.79363,0.0178547,-8.58055e-06,1.67217e-09,-1.1814e-13,-11511.2,-14.0091], Tmin=(864.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.6492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC([O])O(2644)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,D} {8,S}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-55.1457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,1315.93],'cm^-1')),
        HinderedRotor(inertia=(0.140243,'amu*angstrom^2'), symmetry=1, barrier=(3.22447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136697,'amu*angstrom^2'), symmetry=1, barrier=(3.14292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17175,'amu*angstrom^2'), symmetry=1, barrier=(26.9409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74948,0.0552528,-8.37762e-05,7.72178e-08,-2.92091e-11,-6557.02,21.5121], Tmin=(100,'K'), Tmax=(757.128,'K')), NASAPolynomial(coeffs=[5.38536,0.0282609,-1.48808e-05,2.9765e-09,-2.1176e-13,-6884.5,6.45595], Tmin=(757.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.1457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]COC([O])=O(1737)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u1 p0 c0 {4,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {3,D}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-189.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82091,0.0432149,-3.52731e-05,1.43683e-08,-2.34052e-12,-22729,21.3878], Tmin=(100,'K'), Tmax=(1457.83,'K')), NASAPolynomial(coeffs=[12.0197,0.0152316,-6.4807e-06,1.20162e-09,-8.26215e-14,-25702.7,-31.6602], Tmin=(1457.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(CJCO)"""),
)

species(
    label = '[CH]=CO[C](O)O(2645)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {4,S} {8,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u1 p0 c0 {5,D} {9,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-77.2972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,316.503,316.546],'cm^-1')),
        HinderedRotor(inertia=(0.00168226,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194293,'amu*angstrom^2'), symmetry=1, barrier=(13.8159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194337,'amu*angstrom^2'), symmetry=1, barrier=(13.816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194292,'amu*angstrom^2'), symmetry=1, barrier=(13.8156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16452,0.0595778,-7.19529e-05,4.18451e-08,-9.36481e-12,-9192.16,23.1347], Tmin=(100,'K'), Tmax=(1103.79,'K')), NASAPolynomial(coeffs=[14.9903,0.00947502,-3.86588e-06,7.2219e-10,-5.08387e-14,-12244.3,-44.9323], Tmin=(1103.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.2972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = 'C[CH]OC([O])=O(1705)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u1 p0 c0 {1,S} {4,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {3,D}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-207.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3897.47,'J/mol'), sigma=(6.0413,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=608.78 K, Pc=40.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57393,0.0478164,-4.3874e-05,1.99056e-08,-3.5539e-12,-24843.6,21.1131], Tmin=(100,'K'), Tmax=(1355.87,'K')), NASAPolynomial(coeffs=[13.3481,0.0130811,-5.44633e-06,1.01113e-09,-7.00807e-14,-28036.4,-39.275], Tmin=(1355.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(CCsJOC(O))"""),
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
    label = 'C=CO[C]=O(655)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,D} {6,S}
4 C u0 p0 c0 {3,D} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {2,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-97.7509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42599,'amu*angstrom^2'), symmetry=1, barrier=(32.7863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42508,'amu*angstrom^2'), symmetry=1, barrier=(32.7654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06241,0.0351365,-1.90371e-05,-6.74068e-09,6.69218e-12,-11680.2,14.7084], Tmin=(100,'K'), Tmax=(965.056,'K')), NASAPolynomial(coeffs=[13.1388,0.00674637,-2.14126e-06,3.99144e-10,-3.05797e-14,-14633.9,-42.5618], Tmin=(965.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.7509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    label = 'C=COC([O])=O(1617)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {6,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,D} {7,S}
5 C u0 p0 c0 {4,D} {8,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {3,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-287.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,522.833,522.871,522.883,522.888,522.97,523.229,523.239,523.277],'cm^-1')),
        HinderedRotor(inertia=(0.0880078,'amu*angstrom^2'), symmetry=1, barrier=(17.0737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880693,'amu*angstrom^2'), symmetry=1, barrier=(17.0686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74628,0.039296,-1.55324e-05,-1.64506e-08,1.11604e-11,-34528,18.6161], Tmin=(100,'K'), Tmax=(964.027,'K')), NASAPolynomial(coeffs=[15.6906,0.00592718,-1.71663e-06,3.4667e-10,-2.92818e-14,-38354.5,-54.0484], Tmin=(964.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + radical(OC=OOJ)"""),
)

species(
    label = 'C=[C]OC(=O)O(1649)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {6,S}
2 O u0 p2 c0 {4,S} {9,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 C u0 p0 c0 {6,D} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {5,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-291.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6395,0.0435325,-3.02472e-05,1.63766e-09,3.82306e-12,-35007.7,21.653], Tmin=(100,'K'), Tmax=(1015.01,'K')), NASAPolynomial(coeffs=[14.9642,0.00820331,-3.42813e-06,6.99817e-10,-5.36033e-14,-38597.7,-47.1892], Tmin=(1015.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-291.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC(=O)O(2646)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 C u0 p0 c0 {1,S} {6,D} {7,S}
6 C u1 p0 c0 {5,D} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-284.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,375.075,375.082,375.082,375.083,375.089,375.093],'cm^-1')),
        HinderedRotor(inertia=(0.263267,'amu*angstrom^2'), symmetry=1, barrier=(26.2811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263252,'amu*angstrom^2'), symmetry=1, barrier=(26.2814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263236,'amu*angstrom^2'), symmetry=1, barrier=(26.2809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37046,0.0451484,-1.96277e-05,-2.00556e-08,1.41392e-11,-34109.6,19.892], Tmin=(100,'K'), Tmax=(954.205,'K')), NASAPolynomial(coeffs=[19.1316,0.00159789,2.5275e-07,-4.04516e-12,-6.8392e-15,-38906,-72.3349], Tmin=(954.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'O=C(O)O(524)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,423.381,423.383,423.389,423.394],'cm^-1')),
        HinderedRotor(inertia=(0.000940505,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341901,'amu*angstrom^2'), symmetry=1, barrier=(43.4927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98122,0.0114383,2.77968e-05,-4.94729e-08,2.08012e-11,-75136.3,10.9455], Tmin=(100,'K'), Tmax=(956.742,'K')), NASAPolynomial(coeffs=[11.8952,0.00171632,-1.48105e-07,9.26329e-11,-1.38648e-14,-78102.7,-38.2536], Tmin=(956.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs)"""),
)

species(
    label = 'C2H2(23)',
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
    label = 'C[C-]=[O+]C(=O)O(2647)',
    structure = adjacencyList("""1  O u0 p1 c+1 {5,S} {6,D}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u0 p1 c-1 {1,D} {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (125.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1373.69,1562.58,1600,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.130855,'amu*angstrom^2'), symmetry=1, barrier=(3.00861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130855,'amu*angstrom^2'), symmetry=1, barrier=(3.00861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130855,'amu*angstrom^2'), symmetry=1, barrier=(3.00861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97444,0.0572744,-0.000112303,1.18058e-07,-4.52192e-11,15139.2,18.1851], Tmin=(100,'K'), Tmax=(868.843,'K')), NASAPolynomial(coeffs=[0.384686,0.0332469,-1.67033e-05,3.17985e-09,-2.15942e-13,16598.6,32.4399], Tmin=(868.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-O2d)H) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]COC(=O)O(2648)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u0 p1 c0 {4,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-193.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67896,0.0465193,-3.89612e-05,1.56863e-08,-2.50377e-12,-23233.8,20.9299], Tmin=(100,'K'), Tmax=(1487.92,'K')), NASAPolynomial(coeffs=[13.4847,0.014781,-6.96459e-06,1.34991e-09,-9.49198e-14,-26746.9,-40.7173], Tmin=(1487.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-273.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-247.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (53.3569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (74.9161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (153.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (102.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (96.7262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-37.7919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (84.1194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-55.4536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (30.0743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-63.9654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (57.5423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (50.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (46.8649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (54.2169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-153.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (252.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-54.9095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=COC(=O)O(1710)'],
    products = ['CO2(14)', 'CH3CHO(36)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(131.574,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=COC(=O)O(1710)'],
    products = ['O=CCC(=O)O(1703)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(156.914,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O[C]1O[CH]CO1(2640)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.43612e+09,'s^-1'), n=0.0758668, Ea=(56.4791,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=COC([O])[O](2641)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]OC([O])O(2642)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]O[C](O)O(2643)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=COC([O])O(2644)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]COC([O])=O(1737)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CO[C](O)O(2645)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C[CH]OC([O])=O(1705)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.27566e+10,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(=O)O(309)', 'C2H3(29)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C]O(145)', 'CH2CHO(35)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OH(7)', 'C=CO[C]=O(655)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'C=COC([O])=O(1617)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', 'C=[C]OC(=O)O(1649)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.0032e+07,'m^3/(mol*s)'), n=0.204237, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.015586503275756534, var=0.07356775305347633, Tref=1000.0, N=17, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[CH]=COC(=O)O(2646)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.6015e+19,'m^3/(mol*s)'), n=-4.65728, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.23241034472971045, var=0.0, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C
Ea raised from -10.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=COC(=O)O(1710)'],
    products = ['O=C(O)O(524)', 'C2H2(23)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.93334e+10,'s^-1'), n=0.59, Ea=(250.78,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_5CO->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_1BrClFINOPSSi->O_N-3R!H-inRing_5CO->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[C-]=[O+]C(=O)O(2647)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.54807e+17,'s^-1'), n=-1.42952, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5517143818398826, var=16.398565006535087, Tref=1000.0, N=6, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]COC(=O)O(2648)'],
    products = ['C=COC(=O)O(1710)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.40886e+17,'s^-1'), n=-1.43042, Ea=(12.0931,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #577',
    isomers = [
        'C=COC(=O)O(1710)',
    ],
    reactants = [
        ('CO2(14)', 'CH3CHO(36)'),
        ('[O]C(=O)O(309)', 'C2H3(29)'),
        ('O=[C]O(145)', 'CH2CHO(35)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #577',
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

