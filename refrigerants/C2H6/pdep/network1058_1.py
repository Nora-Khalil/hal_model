species(
    label = '[CH2]C(C=C)C(=O)O(1670)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {14,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u0 p0 c0 {4,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-185.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31755,0.0635462,-6.14993e-05,3.52857e-08,-8.81681e-12,-22178.3,24.8523], Tmin=(100,'K'), Tmax=(933.55,'K')), NASAPolynomial(coeffs=[7.84084,0.0355947,-1.6586e-05,3.21096e-09,-2.27022e-13,-23396.3,-6.17002], Tmin=(933.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
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
    label = 'C=C[CH]C(1651)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u1 p0 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800976,'amu*angstrom^2'), symmetry=1, barrier=(1.8416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623781,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]CC=C(54)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u1 p0 c0 {1,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (191.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,345.382],'cm^-1')),
        HinderedRotor(inertia=(0.0750246,'amu*angstrom^2'), symmetry=1, barrier=(6.38,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748783,'amu*angstrom^2'), symmetry=1, barrier=(6.37101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48601,0.0277403,-7.51003e-07,-1.39972e-08,6.14188e-12,23115.6,15.6915], Tmin=(100,'K'), Tmax=(1051,'K')), NASAPolynomial(coeffs=[7.36564,0.0210619,-8.1931e-06,1.49014e-09,-1.03098e-13,21433,-11.2175], Tmin=(1051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""buten3yl1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(O)C=C(2275)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u1 p0 c0 {2,S} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (22.2606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,316.533],'cm^-1')),
        HinderedRotor(inertia=(0.194031,'amu*angstrom^2'), symmetry=1, barrier=(13.7952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00168291,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193905,'amu*angstrom^2'), symmetry=1, barrier=(13.794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79025,0.0401683,-1.22395e-05,-1.42084e-08,8.68906e-12,2764.37,21.9102], Tmin=(100,'K'), Tmax=(989.742,'K')), NASAPolynomial(coeffs=[12.1514,0.0172659,-6.28278e-06,1.14652e-09,-8.14639e-14,-215.843,-32.664], Tmin=(989.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.2606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO)"""),
)

species(
    label = 'C=CC[CH]C(=O)O(2364)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {14,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-188.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3396,0.0526436,-3.34067e-05,1.03146e-08,-1.28263e-12,-22617.2,27.7272], Tmin=(100,'K'), Tmax=(1833.09,'K')), NASAPolynomial(coeffs=[14.5038,0.0239177,-9.90032e-06,1.76558e-09,-1.16692e-13,-27443.4,-43.7598], Tmin=(1833.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C=CCC(=O)O(1678)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {14,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-243.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41908,0.0445247,1.23873e-06,-2.81541e-08,1.17145e-11,-29165.3,24.3353], Tmin=(100,'K'), Tmax=(1111.11,'K')), NASAPolynomial(coeffs=[13.5921,0.027994,-1.32892e-05,2.669e-09,-1.94649e-13,-33555.1,-43.2562], Tmin=(1111.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + radical(Allyl_P)"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787929,2.90875e-06,-1.47516e-09,2.12833e-13,-30281.6,-0.311362], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81028e-07,1.07243e-10,-5.79385e-15,-29958.6,5.9104], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=CC(=C)[C]=O(2365)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,D} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (143.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.07471,'amu*angstrom^2'), symmetry=1, barrier=(24.7096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07027,'amu*angstrom^2'), symmetry=1, barrier=(24.6077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32704,0.0495299,-3.19523e-05,9.86414e-10,4.13383e-12,17372,19.3902], Tmin=(100,'K'), Tmax=(1030.42,'K')), NASAPolynomial(coeffs=[15.5037,0.0128629,-5.31065e-06,1.04709e-09,-7.76085e-14,13475.4,-54.1597], Tmin=(1030.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'C=CO(177)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {7,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u0 p0 c0 {2,D} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H3OH""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C[C]=O(173)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {2,D} {6,S} {7,S}
4 C u1 p0 c0 {1,D} {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (83.3963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.27992,'amu*angstrom^2'), symmetry=1, barrier=(29.4278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0552,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08335,0.0153204,6.5431e-06,-1.77542e-08,7.39398e-12,10067.5,11.895], Tmin=(100,'K'), Tmax=(1002.98,'K')), NASAPolynomial(coeffs=[6.97834,0.0111886,-4.32953e-06,8.0676e-10,-5.75273e-14,8712.69,-9.76666], Tmin=(1002.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CO[C](O)C=C(2366)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {3,S} {14,S}
3  C u1 p0 c0 {1,S} {2,S} {4,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-77.6244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,303.684,303.691,303.692,303.704],'cm^-1')),
        HinderedRotor(inertia=(0.281209,'amu*angstrom^2'), symmetry=1, barrier=(18.4041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28115,'amu*angstrom^2'), symmetry=1, barrier=(18.4037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281203,'amu*angstrom^2'), symmetry=1, barrier=(18.404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281222,'amu*angstrom^2'), symmetry=1, barrier=(18.4036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462037,0.0645886,-3.54806e-05,-1.1734e-08,1.22243e-11,-9196.67,27.2799], Tmin=(100,'K'), Tmax=(956.211,'K')), NASAPolynomial(coeffs=[20.4137,0.0127236,-3.68542e-06,6.55046e-10,-4.95224e-14,-14456.8,-75.6351], Tmin=(956.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.6244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cs_P)"""),
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
    label = 'C=C[CH]C(=O)O(880)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {11,S}
2  O u0 p2 c0 {5,D}
3  C u1 p0 c0 {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {7,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-241.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,546.695,547.306,548.203,548.221,548.301,548.764],'cm^-1')),
        HinderedRotor(inertia=(0.187053,'amu*angstrom^2'), symmetry=1, barrier=(39.6928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186505,'amu*angstrom^2'), symmetry=1, barrier=(39.7058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186312,'amu*angstrom^2'), symmetry=1, barrier=(39.7035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0812,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88004,0.0350486,6.54767e-06,-3.04471e-08,1.23748e-11,-29008.8,17.7914], Tmin=(100,'K'), Tmax=(1094.21,'K')), NASAPolynomial(coeffs=[13.0632,0.0208116,-1.04605e-05,2.16867e-09,-1.61367e-13,-33051.2,-44.4562], Tmin=(1094.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=CCJCO)"""),
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
    label = '[CH]C(C=C)C(=O)O(2367)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u2 p0 c0 {3,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (52.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50187,0.0573497,-4.87178e-05,2.26604e-08,-4.49345e-12,6406.02,24.8309], Tmin=(100,'K'), Tmax=(1160.56,'K')), NASAPolynomial(coeffs=[9.31245,0.0304297,-1.39245e-05,2.67388e-09,-1.88108e-13,4593.08,-14.0137], Tmin=(1160.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'O=C(O)C1[CH]CC1(2368)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {7,S} {14,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {5,S} {13,S}
7  C u0 p0 c0 {1,S} {2,D} {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-164.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99467,0.0279816,4.69606e-05,-7.69272e-08,2.99625e-11,-19668.4,24.3134], Tmin=(100,'K'), Tmax=(991.111,'K')), NASAPolynomial(coeffs=[12.3451,0.0259632,-1.01515e-05,1.95988e-09,-1.4458e-13,-23672.6,-35.3792], Tmin=(991.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsOs) + ring(Cyclobutane) + radical(CCJCC=O)"""),
)

species(
    label = 'C=CC1CO[C]1O(2369)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {14,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-19.0542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64496,0.038439,2.42326e-05,-6.64665e-08,3.17672e-11,-2194.28,24.4936], Tmin=(100,'K'), Tmax=(889.097,'K')), NASAPolynomial(coeffs=[14.3104,0.0180403,-3.07105e-06,2.8436e-10,-1.47853e-14,-5892.35,-43.2523], Tmin=(889.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.0542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C1CC1C(=O)O(2370)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {14,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u1 p0 c0 {4,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-163.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5787,0.0409695,1.10359e-05,-4.33187e-08,1.95569e-11,-19589.6,24.1033], Tmin=(100,'K'), Tmax=(969.668,'K')), NASAPolynomial(coeffs=[12.8485,0.0240321,-8.47767e-06,1.5269e-09,-1.08401e-13,-23164.5,-37.0836], Tmin=(969.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclopropane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1CC1([O])O(2371)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {14,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-38.3247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27396,0.0405749,3.08547e-05,-7.77582e-08,3.49937e-11,-4493.55,24.1776], Tmin=(100,'K'), Tmax=(946.338,'K')), NASAPolynomial(coeffs=[19.1835,0.0134986,-3.29955e-06,5.97111e-10,-4.92984e-14,-10060.5,-72.7413], Tmin=(946.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.3247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH2]C(C=C)=C(O)O(2372)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,S} {14,S}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u1 p0 c0 {3,S} {9,S} {10,S}
7  C u0 p0 c0 {4,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
"""),
    E0 = (-158.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.03466,'amu*angstrom^2'), symmetry=1, barrier=(23.7888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03407,'amu*angstrom^2'), symmetry=1, barrier=(23.7753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03464,'amu*angstrom^2'), symmetry=1, barrier=(23.7885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03405,'amu*angstrom^2'), symmetry=1, barrier=(23.7748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17755,0.088616,-9.62913e-05,4.86846e-08,-9.02414e-12,-18806.9,27.1159], Tmin=(100,'K'), Tmax=(1557.16,'K')), NASAPolynomial(coeffs=[24.9293,0.00407308,1.98678e-06,-6.00146e-10,4.58522e-14,-24818.1,-103.591], Tmin=(1557.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = 'C=CC(=O)O(846)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 C u0 p0 c0 {3,D} {7,S} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-287.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.216788,'amu*angstrom^2'), symmetry=1, barrier=(4.98438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216296,'amu*angstrom^2'), symmetry=1, barrier=(4.97307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38944,0.0170998,-8.39553e-06,1.63166e-09,-9.99489e-14,-34570.6,12.3568], Tmin=(100,'K'), Tmax=(2312.89,'K')), NASAPolynomial(coeffs=[11.8118,0.00529967,-2.53644e-06,4.59873e-10,-2.91762e-14,-39206.4,-36.9378], Tmin=(2312.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'O=[C]O(174)',
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
    label = 'C=CC(=C)C(=O)O(1876)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {3,S} {7,D} {8,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {3,D} {9,S} {10,S}
7  C u0 p0 c0 {4,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-192.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,381.172,381.252,381.261,381.632,381.682,381.842],'cm^-1')),
        HinderedRotor(inertia=(0.00115603,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102548,'amu*angstrom^2'), symmetry=1, barrier=(10.5978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102485,'amu*angstrom^2'), symmetry=1, barrier=(10.6008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922231,0.0619498,-5.7623e-05,2.72634e-08,-5.11201e-12,-23041.4,23.8184], Tmin=(100,'K'), Tmax=(1290.24,'K')), NASAPolynomial(coeffs=[14.83,0.0188329,-7.49656e-06,1.36313e-09,-9.35264e-14,-26630.3,-46.823], Tmin=(1290.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'OH(8)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C([C]=O)C=C(2026)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (239.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,345.144],'cm^-1')),
        HinderedRotor(inertia=(0.115763,'amu*angstrom^2'), symmetry=1, barrier=(9.79813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115904,'amu*angstrom^2'), symmetry=1, barrier=(9.79796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115932,'amu*angstrom^2'), symmetry=1, barrier=(9.79866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69766,0.0540464,-5.77062e-05,3.67601e-08,-9.98028e-12,28849.4,23.3324], Tmin=(100,'K'), Tmax=(875.282,'K')), NASAPolynomial(coeffs=[7.45679,0.0277279,-1.26042e-05,2.40832e-09,-1.6882e-13,27841.2,-3.68495], Tmin=(875.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C(C=C)C([O])=O(1879)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u0 p0 c0 {4,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (40.5328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4196.76,'J/mol'), sigma=(6.40471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=655.52 K, Pc=36.25 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38219,0.064598,-8.86303e-05,8.23345e-08,-3.17795e-11,4962.41,24.7694], Tmin=(100,'K'), Tmax=(782.472,'K')), NASAPolynomial(coeffs=[3.4969,0.0413989,-2.04086e-05,3.97528e-09,-2.78846e-13,5010.72,17.5093], Tmin=(782.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.5328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C=C([O])O(904)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {9,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 C u1 p0 c0 {3,S} {7,S} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
"""),
    E0 = (-66.8439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.08461,'amu*angstrom^2'), symmetry=1, barrier=(24.9374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08269,'amu*angstrom^2'), symmetry=1, barrier=(24.8931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7849,0.0404099,-2.25519e-05,-9.0513e-09,9.11066e-12,-7952.03,16.4535], Tmin=(100,'K'), Tmax=(928.728,'K')), NASAPolynomial(coeffs=[14.8063,0.00565266,-8.58167e-07,1.00566e-10,-8.28922e-15,-11290.4,-50.3559], Tmin=(928.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.8439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(Allyl_P)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56317,0.0223429,1.87066e-05,-3.93097e-08,1.63982e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.265,'K')), NASAPolynomial(coeffs=[9.82996,0.0151966,-5.22271e-06,9.67654e-10,-7.07858e-14,30607.7,-26.985], Tmin=(974.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C(=O)O(2031)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u1 p0 c0 {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {8,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  C u0 p0 c0 {4,D} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-32.5985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62391,0.0565548,-4.92464e-05,2.47374e-08,-5.47967e-12,-3838.72,24.6366], Tmin=(100,'K'), Tmax=(1029.49,'K')), NASAPolynomial(coeffs=[7.53496,0.0335879,-1.57828e-05,3.0674e-09,-2.17336e-13,-5055.79,-4.0526], Tmin=(1029.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.5985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCJ(C)CO) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([C]=C)C(=O)O(2035)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,D} {3,S}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (52.6694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00596,0.0716709,-9.98712e-05,8.4304e-08,-2.93374e-11,6436.9,26.4289], Tmin=(100,'K'), Tmax=(784.631,'K')), NASAPolynomial(coeffs=[7.09913,0.0343082,-1.63998e-05,3.14868e-09,-2.19012e-13,5674.65,-0.253395], Tmin=(784.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.6694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C(=O)O(2037)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u1 p0 c0 {4,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (61.9238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09408,0.0688836,-8.47633e-05,6.17088e-08,-1.88553e-11,7547.83,26.054], Tmin=(100,'K'), Tmax=(787.517,'K')), NASAPolynomial(coeffs=[8.26526,0.0324585,-1.53818e-05,2.97294e-09,-2.08961e-13,6418.38,-6.82964], Tmin=(787.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_P)"""),
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
    label = '[CH2]C=C(O)O(2373)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {9,S}
2  O u0 p2 c0 {4,S} {10,S}
3  C u0 p0 c0 {4,D} {5,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {3,D}
5  C u1 p0 c0 {3,S} {7,S} {8,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-208.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.466095,0.0579037,-6.20701e-05,3.0816e-08,-5.56185e-12,-24909.2,19.7888], Tmin=(100,'K'), Tmax=(1629.41,'K')), NASAPolynomial(coeffs=[17.5629,0.00220168,1.849e-06,-5.08352e-10,3.77846e-14,-28657.8,-65.4466], Tmin=(1629.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[C](C)C(=O)O(1824)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {14,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-243.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45986,0.0517553,-3.11837e-05,8.97388e-09,-1.0368e-12,-29144.8,24.1847], Tmin=(100,'K'), Tmax=(1946.23,'K')), NASAPolynomial(coeffs=[14.595,0.0247599,-1.03782e-05,1.84726e-09,-1.21375e-13,-34257.7,-47.9314], Tmin=(1946.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C=[C]C(C)C(=O)O(2374)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {14,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {7,D} {12,S} {13,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-157.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57482,0.0580851,-5.04277e-05,2.67951e-08,-6.47165e-12,-18900.6,24.0517], Tmin=(100,'K'), Tmax=(939.351,'K')), NASAPolynomial(coeffs=[6.42788,0.0374197,-1.74284e-05,3.37537e-09,-2.38748e-13,-19812.4,0.942031], Tmin=(939.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C(=O)O(2375)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,S} {14,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u1 p0 c0 {5,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
"""),
    E0 = (-148.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50821,0.0575225,-4.49409e-05,1.92598e-08,-3.57152e-12,-17783.6,24.2046], Tmin=(100,'K'), Tmax=(1215.3,'K')), NASAPolynomial(coeffs=[8.85975,0.0333261,-1.50764e-05,2.87735e-09,-2.01508e-13,-19570.5,-12.6959], Tmin=(1215.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)C([O])=O(1821)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u0 p0 c0 {5,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-169.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46597,0.0455907,-2.23275e-05,4.43787e-09,-3.1185e-13,-20400.4,20.4885], Tmin=(100,'K'), Tmax=(2567.09,'K')), NASAPolynomial(coeffs=[28.2222,0.0111726,-5.55562e-06,9.49445e-10,-5.65777e-14,-35507.1,-131.72], Tmin=(2567.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ)"""),
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
    E0 = (-67.5354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (121.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-5.92003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-43.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (124.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (80.1677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (293.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (54.1064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (158.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (283.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-41.5904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (20.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-128.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-19.0075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (68.4943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (20.7085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-42.4027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (42.4925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (286.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (271.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (238.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (101.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (198.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (283.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (293.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (92.5248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-67.3879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (13.3549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-84.9612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-82.1752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['CO2(15)', 'C=C[CH]C(1651)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(98.3198,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(14)', '[CH2]C(O)C=C(2275)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(198.573,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['C=CC[CH]C(=O)O(2364)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['[CH2]C=CCC(=O)O(1678)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO_O] for rate rule [cCs(-HC)CJ;CsJ-HH;CO_O]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CO2(15)', '[CH2]CC=C(54)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(20.4811,'m^3/(mol*s)'), n=1.86833, Ea=(316.938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO2_Cdd;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,3_Insertion_CO2"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H2O(3)', 'C=CC(=C)[C]=O(2365)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.000454348,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_Cdd;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CO(177)', 'C=C[C]=O(173)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;R_OH] for rate rule [Cd_Cdd;Cd_pri_OH]
Euclidian distance = 2.23606797749979
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CO[C](O)C=C(2366)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(112.414,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(T)(17)', 'C=C[CH]C(=O)O(880)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', '[CH]C(C=C)C(=O)O(2367)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['O=C(O)C1[CH]CC1(2368)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['C=CC1CO[C]1O(2369)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['[CH2]C1CC1C(=O)O(2370)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['C=CC1CC1([O])O(2371)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.34239e+09,'s^-1'), n=0.889391, Ea=(146.848,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 143.5 to 146.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=C)=C(O)O(2372)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(410000,'s^-1'), n=2.37, Ea=(207.277,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C2H3(28)', 'C=CC(=O)O(846)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000157915,'m^3/(mol*s)'), n=2.86689, Ea=(2.62026,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.19286017313391585, var=0.167178048265849, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Sp-5R!H-1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Sp-5R!H-1R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]O(174)', 'C=CC=C(61)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(306.062,'m^3/(mol*s)'), n=1.16366, Ea=(33.918,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(6)', 'C=CC(=C)C(=O)O(1876)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1348.96,'m^3/(mol*s)'), n=1.432, Ea=(3.90872,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OH(8)', '[CH2]C([C]=O)C=C(2026)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(6)', '[CH2]C(C=C)C([O])=O(1879)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C2H3(28)', '[CH2]C=C([O])O(904)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=[C]O(174)', '[CH2]C=C[CH2](74)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.63958e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(6)', '[CH2][C](C=C)C(=O)O(2031)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(6)', '[CH2]C([C]=C)C(=O)O(2035)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(6)', '[CH]=CC([CH2])C(=O)O(2037)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.1711e+07,'m^3/(mol*s)'), n=0.21519, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002075092942954368, var=8.990172124599921e-08, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_Sp-3C-2C_N-Sp-4C=3C_Sp-4C-3C"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['C2H2(22)', '[CH2]C=C(O)O(2373)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(23.1678,'s^-1'), n=3.09547, Ea=(258.38,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7BrCClFIOPSSi->O',), comment="""Estimated from node Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7BrCClFIOPSSi->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['C=C[C](C)C(=O)O(1824)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(C)C(=O)O(2374)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=CC(C)C(=O)O(2375)'],
    products = ['[CH2]C(C=C)C(=O)O(1670)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C(=O)O(1670)'],
    products = ['C=CC(C)C([O])=O(1821)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1058',
    isomers = [
        '[CH2]C(C=C)C(=O)O(1670)',
    ],
    reactants = [
        ('CO2(15)', 'C=C[CH]C(1651)'),
        ('CO2(15)', '[CH2]CC=C(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1058',
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

