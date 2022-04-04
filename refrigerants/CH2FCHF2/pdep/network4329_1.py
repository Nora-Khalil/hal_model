species(
    label = 'O=C=C(F)C=C(F)F(7843)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {6,S} {7,D} {9,S}
6 C u0 p0 c0 {1,S} {5,S} {8,D}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 C u0 p0 c0 {4,D} {6,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-514.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,182,240,577,636,1210,1413,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(1.43018,'amu*angstrom^2'), symmetry=1, barrier=(32.8826,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32121,0.0643068,-0.000102582,8.80277e-08,-3.00025e-11,-61766.7,22.3382], Tmin=(100,'K'), Tmax=(797.233,'K')), NASAPolynomial(coeffs=[8.65453,0.0215853,-1.10486e-05,2.15879e-09,-1.50723e-13,-62747.6,-10.1976], Tmin=(797.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-514.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cd-Cd(CCO)H) + group(Cd(Cdd-Od)CF) + group(CdCFF) + missing(Cdd-CdO2d)"""),
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
    label = 'F[C]C=C(F)F-2(6976)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {5,D} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 C u0 p1 c0 {3,S} {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-223.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,315,622,1128],'cm^-1')),
        HinderedRotor(inertia=(0.0539877,'amu*angstrom^2'), symmetry=1, barrier=(22.3783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53736,0.0330301,-3.01871e-05,1.357e-08,-2.48101e-12,-26862.3,15.9167], Tmin=(100,'K'), Tmax=(1281.09,'K')), NASAPolynomial(coeffs=[8.95577,0.0129898,-6.72251e-06,1.35932e-09,-9.81616e-14,-28506.9,-16.6384], Tmin=(1281.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCsH) + group(CdCFF) + group(CJ2_singlet-FC)"""),
)

species(
    label = 'O=C=[C]F(289)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (33.6712,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(58.9933,'amu')),
        NonlinearRotor(inertia=([3.58812,117.154,120.742],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([280.764,365.894,563.98,866.433,1429.66,2068.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0191,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94745,0.00353157,4.23087e-05,-1.10022e-07,8.19717e-11,4052.13,8.19715], Tmin=(10,'K'), Tmax=(472.411,'K')), NASAPolynomial(coeffs=[4.41104,0.00925722,-6.51505e-06,2.12227e-09,-2.60232e-13,3900.64,5.16846], Tmin=(472.411,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(33.6712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF2CH(73)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u1 p0 c0 {3,D} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-92.0165,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(63.0046,'amu')),
        NonlinearRotor(inertia=([43.4896,46.129,89.6186],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([420.994,524.562,552.014,640.304,746.121,978.508,1261.79,1779.04,3347.18],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0261,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95511,0.00251604,4.75146e-05,-8.97685e-08,4.98372e-11,-11065.2,8.58309], Tmin=(10,'K'), Tmax=(611.191,'K')), NASAPolynomial(coeffs=[3.6485,0.0155109,-1.13452e-05,3.84899e-09,-4.87703e-13,-11233,8.2324], Tmin=(611.191,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-92.0165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[CH]DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=C=CC(F)(F)F(15435)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {9,S}
7 C u0 p0 c0 {6,D} {8,D}
8 C u0 p0 c0 {4,D} {7,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-592.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71227,0.0538806,-6.50211e-05,4.11727e-08,-1.05887e-11,-71128.4,20.2914], Tmin=(100,'K'), Tmax=(937.906,'K')), NASAPolynomial(coeffs=[9.794,0.0194123,-9.89368e-06,1.98657e-09,-1.43189e-13,-72644.3,-18.1799], Tmin=(937.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-592.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCdFFF) + group(Cds-CdsCsH) + group(Cdd-(Cdd-O2d)Cds) + missing(Cdd-CddO2d)"""),
)

species(
    label = 'O=C1C(F)=CC1(F)F(15436)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {5,S} {8,D} {9,S}
7 C u0 p0 c0 {4,D} {5,S} {8,S}
8 C u0 p0 c0 {3,S} {6,D} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-572.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07265,0.048862,-5.53336e-05,2.9394e-08,-1.43011e-12,-68806.3,18.4867], Tmin=(100,'K'), Tmax=(534.991,'K')), NASAPolynomial(coeffs=[5.79702,0.0277934,-1.52648e-05,3.14334e-09,-2.28981e-13,-69301.8,1.9419], Tmin=(534.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_cyclic(Cs(F)2-CO) + longDistanceInteraction_cyclic(Cd(F)-CO) + ring(Cyclobutene)"""),
)

species(
    label = 'O=[C]C(F)[C]=C(F)F(10514)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u1 p0 c0 {5,S} {6,D}
8 C u1 p0 c0 {4,D} {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-281.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,562,600,623,1070,1265,1685,370,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.711412,'amu*angstrom^2'), symmetry=1, barrier=(16.3568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710501,'amu*angstrom^2'), symmetry=1, barrier=(16.3358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.101,0.0721779,-0.000129327,1.18221e-07,-4.13953e-11,-33724.4,25.4812], Tmin=(100,'K'), Tmax=(842.363,'K')), NASAPolynomial(coeffs=[8.07503,0.0226949,-1.2068e-05,2.35427e-09,-1.62555e-13,-34318.7,-3.52153], Tmin=(842.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFF) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC(F)=[C][C](F)F(10281)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {6,S} {8,D}
6 C u0 p0 c0 {4,D} {5,S} {9,S}
7 C u1 p0 c0 {2,S} {3,S} {8,S}
8 C u1 p0 c0 {5,D} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-312.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,384,691,1241,2782.5,750,1395,475,1775,1000,161,297,490,584,780,1358,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.197442,'amu*angstrom^2'), symmetry=1, barrier=(4.53958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0277355,'amu*angstrom^2'), symmetry=1, barrier=(33.7222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50096,0.0610025,-9.06863e-05,7.62157e-08,-2.66034e-11,-37456,23.2286], Tmin=(100,'K'), Tmax=(694.019,'K')), NASAPolynomial(coeffs=[7.61835,0.0257382,-1.44546e-05,2.97498e-09,-2.15686e-13,-38305,-4.04879], Tmin=(694.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(Cds-CdsCsH) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)H) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cds_S)"""),
)

species(
    label = 'O[C]=C(F)[C]=C(F)F(15437)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {8,S} {9,S}
5 C u0 p0 c0 {1,S} {7,S} {8,D}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u1 p0 c0 {5,S} {6,D}
8 C u1 p0 c0 {4,S} {5,D}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-209.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,250,446,589,854,899,562,600,623,1070,1265,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.580551,'amu*angstrom^2'), symmetry=1, barrier=(15.0682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15891,'amu*angstrom^2'), symmetry=1, barrier=(38.1618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22615,0.0602378,-7.11959e-05,4.07195e-08,-9.08239e-12,-25108.5,27.3782], Tmin=(100,'K'), Tmax=(1097.71,'K')), NASAPolynomial(coeffs=[14.2839,0.0126559,-6.17579e-06,1.23093e-09,-8.89431e-14,-27975.2,-36.835], Tmin=(1097.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(CdCFF) + radical(Cdj(Cd-F1sCd)(Cd-F1sF1s)) + radical(C=CJO)"""),
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
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.158,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C=C=C[C](F)F(15438)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {5,S} {6,D} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u0 p0 c0 {4,D} {7,D}
7 C u0 p0 c0 {3,D} {6,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-172.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,161,297,490,584,780,1358,540,610,2055,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(1.22486,'amu*angstrom^2'), symmetry=1, barrier=(28.162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9829,0.0485414,-6.5204e-05,4.80577e-08,-1.4553e-11,-20698.3,20.5725], Tmin=(100,'K'), Tmax=(798.834,'K')), NASAPolynomial(coeffs=[7.84856,0.0191686,-1.00463e-05,2.02312e-09,-1.45401e-13,-21635.4,-6.40808], Tmin=(798.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cdd-(Cdd-O2d)Cds) + missing(Cdd-CddO2d) + radical(CsCdF1sF1s)"""),
)

species(
    label = 'O=C=C(F)C=[C]F(11077)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {5,S} {6,D} {8,S}
5 C u0 p0 c0 {1,S} {4,S} {7,D}
6 C u1 p0 c0 {2,S} {4,D}
7 C u0 p0 c0 {3,D} {5,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-57.5744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,167,640,1190,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(1.27041,'amu*angstrom^2'), symmetry=1, barrier=(29.2091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5152,0.063353,-0.000117941,1.10351e-07,-3.88753e-11,-6843.46,20.563], Tmin=(100,'K'), Tmax=(864.244,'K')), NASAPolynomial(coeffs=[6.58686,0.0207439,-1.0775e-05,2.06526e-09,-1.40424e-13,-7005.44,0.969661], Tmin=(864.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.5744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cd-Cd(CCO)H) + group(Cd(Cdd-Od)CF) + group(CdCFH) + missing(Cdd-CdO2d) + radical(Cdj(Cd-CdH)(F1s))"""),
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
    label = 'O=[C]C(F)=C=C(F)F(15439)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {7,D} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u0 p0 c0 {5,D} {6,D}
8 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-340.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.43783,'amu*angstrom^2'), symmetry=1, barrier=(33.0586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59031,0.0583531,-9.45247e-05,8.32186e-08,-2.95295e-11,-40866.6,22.2968], Tmin=(100,'K'), Tmax=(744.434,'K')), NASAPolynomial(coeffs=[7.90592,0.0207881,-1.15189e-05,2.33396e-09,-1.66768e-13,-41706.4,-5.6329], Tmin=(744.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCddCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)H) + group(CdCddFF) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C=C(F)[C]C(F)F(15440)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
6 C u0 p0 c0 {3,S} {7,S} {8,D}
7 C u0 p1 c0 {5,S} {6,S}
8 C u0 p0 c0 {4,D} {6,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-213.761,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([257,409,1143,1361,2944,3010,987.5,1337.5,450,1655,2120,512.5,787.5,226.212,226.212,226.213,226.214,226.214,226.215],'cm^-1')),
        HinderedRotor(inertia=(1.33443,'amu*angstrom^2'), symmetry=1, barrier=(48.4571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33443,'amu*angstrom^2'), symmetry=1, barrier=(48.4571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.965021,0.0758185,-0.000134386,1.24083e-07,-4.4116e-11,-25608.9,34.1437], Tmin=(100,'K'), Tmax=(833.994,'K')), NASAPolynomial(coeffs=[7.58173,0.0265384,-1.41959e-05,2.78189e-09,-1.92993e-13,-26102.4,7.08115], Tmin=(833.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCFFH) + group(Cd(Cdd-Od)CF) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C(F)[C]C=C(F)F(15441)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {6,D} {8,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 C u0 p0 c0 {3,S} {4,D} {8,S}
8 C u0 p1 c0 {5,S} {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-308.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,525,764,1175,1295,1876,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14338,'amu*angstrom^2'), symmetry=1, barrier=(49.2806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14284,'amu*angstrom^2'), symmetry=1, barrier=(49.268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15502,0.0697931,-0.000116482,1.0484e-07,-3.70402e-11,-36985.7,33.5619], Tmin=(100,'K'), Tmax=(815.315,'K')), NASAPolynomial(coeffs=[7.6815,0.025804,-1.35299e-05,2.65162e-09,-1.84854e-13,-37652,5.84776], Tmin=(815.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCsH) + group(CdCFF) + group(COCFO) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=C=C(F)C(F)[C]F(15442)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6 C u0 p0 c0 {2,S} {5,S} {8,D}
7 C u0 p1 c0 {3,S} {5,S}
8 C u0 p0 c0 {4,D} {6,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-259.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,145,326,398,834,1303,617,898,1187,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(1.01529,'amu*angstrom^2'), symmetry=1, barrier=(23.3436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901318,'amu*angstrom^2'), symmetry=1, barrier=(20.7231,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19063,0.0671389,-0.000107781,9.13537e-08,-3.07566e-11,-31112,24.219], Tmin=(100,'K'), Tmax=(782.002,'K')), NASAPolynomial(coeffs=[9.52668,0.020522,-1.07332e-05,2.11517e-09,-1.48426e-13,-32294.1,-13.1702], Tmin=(782.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cd(Cdd-Od)CF) + group(CJ2_singlet-FCs) + missing(Cdd-CdO2d)"""),
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
    label = 'O=C=C=C=C(F)F(15443)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {4,D} {6,D}
6 C u0 p0 c0 {5,D} {7,D}
7 C u0 p0 c0 {3,D} {6,D}
"""),
    E0 = (-141.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([94,120,354,641,825,1294,540,563.333,586.667,610,1970,2140,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09941,0.0447337,-6.01696e-05,4.18696e-08,-1.1656e-11,-16949.2,17.5441], Tmin=(100,'K'), Tmax=(875.463,'K')), NASAPolynomial(coeffs=[8.88696,0.0137209,-7.03214e-06,1.40477e-09,-1.00625e-13,-18137.6,-14.2991], Tmin=(875.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(CdCddFF) + group(Cdd-CdsCds) + group(Cdd-(Cdd-O2d)Cds) + missing(Cdd-CddO2d)"""),
)

species(
    label = 'O=C=C(F)C#CF(11162)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,D}
5 C u0 p0 c0 {4,S} {7,T}
6 C u0 p0 c0 {3,D} {4,D}
7 C u0 p0 c0 {2,S} {5,T}
"""),
    E0 = (-52.3256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2175,525,2120,512.5,787.5,239,401,1367,180],'cm^-1')),
        HinderedRotor(inertia=(1.58656,'amu*angstrom^2'), symmetry=1, barrier=(36.4782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.955,0.0496171,-8.14164e-05,7.22065e-08,-2.57394e-11,-6224.13,18.942], Tmin=(100,'K'), Tmax=(747.835,'K')), NASAPolynomial(coeffs=[7.306,0.0174455,-9.76542e-06,1.98392e-09,-1.41871e-13,-6925.18,-4.65478], Tmin=(747.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.3256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cd(Cdd-Od)CF) + group(Ct-Ct(Cds-Cds)) + missing(Cdd-CdO2d) + group(CtCF)"""),
)

species(
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = '[CH]C(F)=C=O(15444)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p1 c0 {3,S} {6,S}
5 C u0 p0 c0 {2,D} {3,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (226.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0314657,'amu*angstrom^2'), symmetry=1, barrier=(8.8627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0377,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46209,0.0385916,-6.88992e-05,6.38797e-08,-2.26961e-11,27283,14.6071], Tmin=(100,'K'), Tmax=(840.564,'K')), NASAPolynomial(coeffs=[5.72802,0.0135133,-7.12819e-06,1.39039e-09,-9.61292e-14,27070.8,1.42218], Tmin=(840.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cd(Cdd-Od)CF) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
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
    E0 = (-156.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-61.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-209.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-76.1947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-41.0148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (11.7694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (82.3687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (197.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (123.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (54.7519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-3.27256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (3.88095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (35.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-19.4257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (129.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (204.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C=C(F)C=C(F)F(7843)'],
    products = ['CO(13)', 'F[C]C=C(F)F-2(6976)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(7.15699e+14,'s^-1'), n=0.0573689, Ea=(175.793,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='RFC=C=O',), comment="""Estimated from node RFC=C=O"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C=C(F)C=C(F)F(7843)'],
    products = ['O=C=C=CC(F)(F)F(15435)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.72966e+11,'s^-1'), n=0.63878, Ea=(270.836,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C=C(F)C=C(F)F(7843)'],
    products = ['O=C1C(F)=CC1(F)F(15436)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]C(F)[C]=C(F)F(10514)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=CC(F)=[C][C](F)F(10281)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O[C]=C(F)[C]=C(F)F(15437)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'O=C=C=C[C](F)F(15438)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'O=C=C(F)C=[C]F(11077)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C=[C]F(289)', 'CF2CH(73)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -31.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'O=[C]C(F)=C=C(F)F(15439)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(1.2655,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C=C(F)[C]C(F)F(15440)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.33333e+12,'s^-1'), n=8.2394e-08, Ea=(28.3428,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C(F)[C]C=C(F)F(15441)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(130.046,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C=C(F)C(F)[C]F(15442)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(113.249,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'O=C=C=C=C(F)F(15443)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(221.014,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'O=C=C(F)C#CF(11162)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(281.042,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2(43)', '[CH]C(F)=C=O(15444)'],
    products = ['O=C=C(F)C=C(F)F(7843)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.66533e+06,'cm^3/(mol*s)'), n=1.53, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for CF2
Exact match found for rate rule [CF2]
Euclidian distance = 0
family: halocarbene_recombination_double"""),
)

network(
    label = 'PDepNetwork #4329',
    isomers = [
        'O=C=C(F)C=C(F)F(7843)',
    ],
    reactants = [
        ('CO(13)', 'F[C]C=C(F)F-2(6976)'),
        ('O=C=[C]F(289)', 'CF2CH(73)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4329',
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

