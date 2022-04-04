species(
    label = 'C[CH]OO[CH]C(F)F(6716)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
6  C u0 p0 c0 {8,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {5,S} {13,S}
8  C u1 p0 c0 {4,S} {6,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-260.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,522,611,926,1093,1137,1374,1416,3112,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.111003,'amu*angstrom^2'), symmetry=1, barrier=(2.55218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307254,'amu*angstrom^2'), symmetry=1, barrier=(7.06437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.780751,'amu*angstrom^2'), symmetry=1, barrier=(49.2077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13995,'amu*angstrom^2'), symmetry=1, barrier=(49.2017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13973,'amu*angstrom^2'), symmetry=1, barrier=(49.1967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.991005,0.0784218,-7.09949e-05,-1.46583e-08,5.26342e-11,-31192.8,25.9444], Tmin=(100,'K'), Tmax=(497.11,'K')), NASAPolynomial(coeffs=[8.39338,0.0402498,-2.03601e-05,3.99898e-09,-2.81699e-13,-32193,-7.25242], Tmin=(497.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CCsJOOC) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CC(F)F(990)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-557.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([195,270,1147,1130,1359,1388,1409,3075,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.595284,'amu*angstrom^2'), symmetry=1, barrier=(13.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93683,0.00748989,0.000222302,-1.39306e-06,2.75418e-09,-67109.3,9.32225], Tmin=(10,'K'), Tmax=(169.939,'K')), NASAPolynomial(coeffs=[3.97443,0.0203448,-1.24424e-05,3.60772e-09,-4.02215e-13,-67130.4,8.62375], Tmin=(169.939,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-557.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODCC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.65,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]C(F)F(506)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u2 p0 c0 {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-68.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,200.055,1350.79,1616.28],'cm^-1')),
        HinderedRotor(inertia=(0.00421209,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96601,0.019505,-1.18279e-05,2.09012e-09,3.19809e-13,-8190.14,12.8831], Tmin=(100,'K'), Tmax=(1211.05,'K')), NASAPolynomial(coeffs=[7.87164,0.00800482,-3.40878e-06,6.62037e-10,-4.73286e-14,-9723.2,-13.1469], Tmin=(1211.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]O[O](416)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (163.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.385356,'amu*angstrom^2'), symmetry=1, barrier=(8.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119494,'amu*angstrom^2'), symmetry=1, barrier=(8.86708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06788,0.0231601,-1.22439e-05,-6.93875e-09,1.03576e-11,19636.5,14.5219], Tmin=(100,'K'), Tmax=(545.313,'K')), NASAPolynomial(coeffs=[4.49937,0.0169793,-7.12433e-06,1.32857e-09,-9.2142e-14,19416.2,7.89488], Tmin=(545.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[O]O[CH]C(F)F(786)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-238.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.648643,'amu*angstrom^2'), symmetry=1, barrier=(14.9136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647571,'amu*angstrom^2'), symmetry=1, barrier=(14.8889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39956,0.0613837,-0.000102358,8.44858e-08,-2.68342e-11,-28585.3,19.1592], Tmin=(100,'K'), Tmax=(853.674,'K')), NASAPolynomial(coeffs=[10.6471,0.0122374,-5.78324e-06,1.08643e-09,-7.3458e-14,-29952.2,-22.7506], Tmin=(853.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH]C-2(510)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u2 p0 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82364,-0.000909745,3.21388e-05,-3.73491e-08,1.33095e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.803,'K')), NASAPolynomial(coeffs=[4.3048,0.0094308,-3.27565e-06,5.95137e-10,-4.2732e-14,40709.2,1.84239], Tmin=(960.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC1OOC1C(F)F(6717)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
8  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-498.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16009,0.063248,-5.09153e-05,2.13069e-08,-3.71889e-12,-59902.3,19.4091], Tmin=(100,'K'), Tmax=(1316.53,'K')), NASAPolynomial(coeffs=[11.7569,0.0310525,-1.42337e-05,2.73232e-09,-1.91767e-13,-62692.6,-34.6285], Tmin=(1316.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-498.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsFFH) + group(Cs-CsHHH) + ring(12dioxetane)"""),
)

species(
    label = 'CCOOC=C(F)F(6720)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {8,D} {14,S}
8  C u0 p0 c0 {1,S} {2,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-437.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06873,0.0709189,-7.4508e-05,4.61117e-08,-1.24167e-11,-52479.2,26.4315], Tmin=(100,'K'), Tmax=(870.415,'K')), NASAPolynomial(coeffs=[8.03772,0.0388928,-1.93166e-05,3.8395e-09,-2.75254e-13,-53692.4,-6.22264], Tmin=(870.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-437.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'C=COOCC(F)F(6718)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
7  C u0 p0 c0 {4,S} {8,D} {12,S}
8  C u0 p0 c0 {7,D} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-490.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671966,0.0745479,-7.36606e-05,3.77498e-08,-7.88142e-12,-58828.3,27.0234], Tmin=(100,'K'), Tmax=(1140.18,'K')), NASAPolynomial(coeffs=[13.6696,0.028949,-1.36709e-05,2.6734e-09,-1.90366e-13,-61792.2,-37.3877], Tmin=(1140.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-490.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[CH][O](3430)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123964,'amu*angstrom^2'), symmetry=1, barrier=(2.85017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65563,0.0114443,2.34958e-06,-4.83183e-09,1.17971e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.62,'K')), NASAPolynomial(coeffs=[6.06246,0.0136329,-6.35985e-06,1.18413e-09,-7.90693e-14,16986.2,-5.89951], Tmin=(1718.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C(F)F(1636)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-255.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.637092,'amu*angstrom^2'), symmetry=1, barrier=(14.648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0921,0.0467254,-7.94794e-05,7.13576e-08,-2.49602e-11,-30696.7,15.5487], Tmin=(100,'K'), Tmax=(816.59,'K')), NASAPolynomial(coeffs=[6.94864,0.0153515,-7.91621e-06,1.55898e-09,-1.09043e-13,-31237,-5.34886], Tmin=(816.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'F(37)',
    structure = adjacencyList("""multiplicity 2
1 F u1 p3 c0
"""),
    E0 = (72.8916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.9984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C[CH]OOC=CF(6818)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {4,S} {11,S}
6  C u0 p0 c0 {3,S} {7,D} {12,S}
7  C u0 p0 c0 {1,S} {6,D} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-58.3287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(1.73895,'amu*angstrom^2'), symmetry=1, barrier=(39.9819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47891,'amu*angstrom^2'), symmetry=1, barrier=(11.0111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00620655,'amu*angstrom^2'), symmetry=1, barrier=(39.9816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73891,'amu*angstrom^2'), symmetry=1, barrier=(39.9809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41603,0.0625833,-6.34838e-05,3.84857e-08,-1.02786e-11,-6927.19,24.4005], Tmin=(100,'K'), Tmax=(871.828,'K')), NASAPolynomial(coeffs=[7.20841,0.0360072,-1.77585e-05,3.52034e-09,-2.52052e-13,-7937.18,-2.74988], Tmin=(871.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.3287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CCsJOOC)"""),
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
    label = 'C[CH]OOC=C(F)F(6819)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {5,S} {12,S}
7  C u0 p0 c0 {4,S} {8,D} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-250.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(1.82733,'amu*angstrom^2'), symmetry=1, barrier=(42.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333024,'amu*angstrom^2'), symmetry=1, barrier=(7.65688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00804854,'amu*angstrom^2'), symmetry=1, barrier=(42.0419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82649,'amu*angstrom^2'), symmetry=1, barrier=(41.9945,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0678,0.0724024,-8.43123e-05,5.08152e-08,-8.70072e-12,-30057.5,26.7537], Tmin=(100,'K'), Tmax=(585.549,'K')), NASAPolynomial(coeffs=[7.53702,0.0375654,-1.90364e-05,3.78256e-09,-2.69899e-13,-30975.5,-2.36376], Tmin=(585.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CCsJOOC)"""),
)

species(
    label = 'C=COO[CH]C(F)F(6820)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u1 p0 c0 {3,S} {5,S} {10,S}
7  C u0 p0 c0 {4,S} {8,D} {11,S}
8  C u0 p0 c0 {7,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-303.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.45114,'amu*angstrom^2'), symmetry=1, barrier=(33.3647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766083,'amu*angstrom^2'), symmetry=1, barrier=(17.6138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765876,'amu*angstrom^2'), symmetry=1, barrier=(17.609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94334,'amu*angstrom^2'), symmetry=1, barrier=(44.6812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751197,0.0754287,-8.39322e-05,4.94113e-08,-1.18843e-11,-36409.6,27.043], Tmin=(100,'K'), Tmax=(997.028,'K')), NASAPolynomial(coeffs=[12.2953,0.0291145,-1.42539e-05,2.82064e-09,-2.01957e-13,-38711.6,-28.6164], Tmin=(997.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOOC)"""),
)

species(
    label = 'C[CH]OOC[C](F)F(6821)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
7  C u1 p0 c0 {4,S} {6,S} {14,S}
8  C u1 p0 c0 {1,S} {2,S} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-244.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,198.088,1796.86],'cm^-1')),
        HinderedRotor(inertia=(0.381207,'amu*angstrom^2'), symmetry=1, barrier=(9.59733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35017,'amu*angstrom^2'), symmetry=1, barrier=(36.9939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378426,'amu*angstrom^2'), symmetry=1, barrier=(9.62405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38815,'amu*angstrom^2'), symmetry=1, barrier=(36.9917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.34574,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.18121,0.0921779,-0.000139335,1.2128e-07,-4.28433e-11,-29231.6,29.2303], Tmin=(100,'K'), Tmax=(777.501,'K')), NASAPolynomial(coeffs=[8.67474,0.0393982,-1.99858e-05,3.91817e-09,-2.75156e-13,-30277.8,-7.84294], Tmin=(777.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(CsCsFFH) + radical(CCsJOOC) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]COO[CH]C(F)F(6822)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u1 p0 c0 {4,S} {6,S} {12,S}
8  C u1 p0 c0 {5,S} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-232.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.204834,0.0906691,-0.000129936,1.0676e-07,-3.61133e-11,-27848.7,30.3091], Tmin=(100,'K'), Tmax=(729.174,'K')), NASAPolynomial(coeffs=[9.75024,0.0374976,-1.88921e-05,3.71413e-09,-2.6217e-13,-29219.3,-12.5798], Tmin=(729.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CCsJOOC) + radical(CJCOOH)"""),
)

species(
    label = 'CCOO[CH][C](F)F(6823)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  C u1 p0 c0 {4,S} {8,S} {14,S}
8  C u1 p0 c0 {1,S} {2,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-244.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,198.088,1796.86],'cm^-1')),
        HinderedRotor(inertia=(0.381207,'amu*angstrom^2'), symmetry=1, barrier=(9.59733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35017,'amu*angstrom^2'), symmetry=1, barrier=(36.9939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378426,'amu*angstrom^2'), symmetry=1, barrier=(9.62405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38815,'amu*angstrom^2'), symmetry=1, barrier=(36.9917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.34574,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.18121,0.0921779,-0.000139335,1.2128e-07,-4.28433e-11,-29231.6,29.2303], Tmin=(100,'K'), Tmax=(777.501,'K')), NASAPolynomial(coeffs=[8.67474,0.0393982,-1.99858e-05,3.91817e-09,-2.75156e-13,-30277.8,-7.84294], Tmin=(777.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(CsCsFFH) + radical(CCsJOOC) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2][CH]OOCC(F)F(6824)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
7  C u1 p0 c0 {4,S} {8,S} {12,S}
8  C u1 p0 c0 {7,S} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-232.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,258.082,1749.3],'cm^-1')),
        HinderedRotor(inertia=(0.176738,'amu*angstrom^2'), symmetry=1, barrier=(8.35359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176739,'amu*angstrom^2'), symmetry=1, barrier=(8.35359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863523,'amu*angstrom^2'), symmetry=1, barrier=(40.8146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863524,'amu*angstrom^2'), symmetry=1, barrier=(40.8145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863524,'amu*angstrom^2'), symmetry=1, barrier=(40.8145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.204834,0.0906691,-0.000129936,1.0676e-07,-3.61133e-11,-27848.7,30.3091], Tmin=(100,'K'), Tmax=(729.174,'K')), NASAPolynomial(coeffs=[9.75024,0.0374976,-1.88921e-05,3.71413e-09,-2.6217e-13,-29219.3,-12.5798], Tmin=(729.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cs-CsHHH) + radical(CCsJOOC) + radical(CJCOOH)"""),
)

species(
    label = 'C[CH]OOC(F)[CH]F(6825)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {6,S} {13,S}
8  C u1 p0 c0 {2,S} {5,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-237.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,487,638,688,1119,1325,1387,3149,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,334,575,1197,1424,3202,193.811,193.811],'cm^-1')),
        HinderedRotor(inertia=(0.176587,'amu*angstrom^2'), symmetry=1, barrier=(4.70712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58094,'amu*angstrom^2'), symmetry=1, barrier=(42.1418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00351,'amu*angstrom^2'), symmetry=1, barrier=(26.7487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58101,'amu*angstrom^2'), symmetry=1, barrier=(42.1418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.48802,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394807,0.0864402,-0.00010941,7.78314e-08,-2.30032e-11,-28438.8,27.7538], Tmin=(100,'K'), Tmax=(814.447,'K')), NASAPolynomial(coeffs=[10.3981,0.0373115,-1.89288e-05,3.76909e-09,-2.69562e-13,-30068.3,-18.4533], Tmin=(814.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CCsJOOC) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'CC(F)OO[CH][CH]F(6826)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {8,S} {13,S}
8  C u1 p0 c0 {2,S} {7,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-256.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,261,493,600,1152,1365,1422,3097,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,334,575,1197,1424,3202,180,2433.69],'cm^-1')),
        HinderedRotor(inertia=(1.54068,'amu*angstrom^2'), symmetry=1, barrier=(35.4232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54739,'amu*angstrom^2'), symmetry=1, barrier=(35.5776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55073,'amu*angstrom^2'), symmetry=1, barrier=(35.6543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53909,'amu*angstrom^2'), symmetry=1, barrier=(35.3866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580153,'amu*angstrom^2'), symmetry=1, barrier=(13.3389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.086,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270004,0.0897194,-0.000124217,9.89191e-08,-3.27242e-11,-30770,28.4397], Tmin=(100,'K'), Tmax=(731.593,'K')), NASAPolynomial(coeffs=[9.55996,0.0389224,-2.00591e-05,3.99734e-09,-2.84998e-13,-32129.2,-13.4749], Tmin=(731.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(CsCsFHH) + radical(CCsJOOC) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
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
    E0 = (-110.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (244.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (255.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-102.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-102.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-102.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-100.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-107.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (204.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (110.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (61.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (59.6636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (71.6913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (77.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-58.8078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-36.4547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (89.3541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (89.3001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]OO[CH]C(F)F(6716)'],
    products = ['O=CC(F)F(990)', 'CH3CHO(36)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(F)F(506)', 'C[CH]O[O](416)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]O[CH]C(F)F(786)', '[CH]C-2(510)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C[CH]OO[CH]C(F)F(6716)'],
    products = ['CC1OOC1C(F)F(6717)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C[CH]OO[CH]C(F)F(6716)'],
    products = ['CCOOC=C(F)F(6720)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C[CH]OO[CH]C(F)F(6716)'],
    products = ['C=COOCC(F)F(6718)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=CC(F)F(990)', 'C[CH][O](3430)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(149.753,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][CH]C(F)F(1636)', 'CH3CHO(36)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(177.521,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'C[CH]OOC=CF(6818)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.8984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'C[CH]OOC=C(F)F(6819)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.18599,'m^3/(mol*s)'), n=2.09886, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.020567443735397595, var=1.0626053141426195, Tref=1000.0, N=111, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C=COO[CH]C(F)F(6820)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.67e+06,'m^3/(mol*s)'), n=0.1, Ea=(3.42452,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_4R!H->O_Sp-4O-2C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_Sp-2CS=1CCSS_1CS->C_2CS->C_Ext-2C-R_4R!H->O_Sp-4O-2C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][CH]C(F)F(1636)', 'C[CH][O](3430)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(8.09286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]OOC[C](F)F(6821)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.82764e+09,'s^-1'), n=1.1605, Ea=(166.084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out] for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]OO[CH]C(F)F(6716)'],
    products = ['[CH2]COO[CH]C(F)F(6822)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.41e+13,'s^-1','+|-',2), n=0, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1500,'K'), comment="""From training reaction 346 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCOO[CH][C](F)F(6823)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(121.56,'s^-1'), n=2.45, Ea=(35.5849,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_noH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]OOCC(F)F(6824)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]OOC(F)[CH]F(6825)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(177.104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC(F)OO[CH][CH]F(6826)'],
    products = ['C[CH]OO[CH]C(F)F(6716)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(196.463,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #2498',
    isomers = [
        'C[CH]OO[CH]C(F)F(6716)',
    ],
    reactants = [
        ('O=CC(F)F(990)', 'CH3CHO(36)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2498',
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

