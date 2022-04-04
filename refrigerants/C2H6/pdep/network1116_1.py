species(
    label = 'C=C[CH]CO[C]=O(1880)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u1 p0 c0 {3,S} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-12.7977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,649.33,649.632,649.9],'cm^-1')),
        HinderedRotor(inertia=(1.87626,'amu*angstrom^2'), symmetry=1, barrier=(43.1389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763492,'amu*angstrom^2'), symmetry=1, barrier=(22.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990912,'amu*angstrom^2'), symmetry=1, barrier=(22.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.827254,'amu*angstrom^2'), symmetry=1, barrier=(22.7822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815359,0.0612583,-4.82232e-05,1.68331e-08,-1.73287e-12,-1417.17,24.7599], Tmin=(100,'K'), Tmax=(1151.52,'K')), NASAPolynomial(coeffs=[15.5676,0.0203457,-8.38789e-06,1.56239e-09,-1.09173e-13,-5499.67,-51.4669], Tmin=(1151.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.7977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=CCJCO) + radical((O)CJOCC)"""),
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
    label = '[C]=O(152)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08918,0.00200392,-1.61651e-05,2.55044e-08,-1.16417e-11,52802.7,4.52499], Tmin=(100,'K'), Tmax=(856.118,'K')), NASAPolynomial(coeffs=[0.961586,0.00569052,-3.48048e-06,7.19212e-10,-5.0805e-14,53738.7,21.4665], Tmin=(856.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C=C[CH]C[O](1223)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3  C u1 p0 c0 {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (166.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,308.786,310.533,1588.88],'cm^-1')),
        HinderedRotor(inertia=(0.595621,'amu*angstrom^2'), symmetry=1, barrier=(40.4269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595433,'amu*angstrom^2'), symmetry=1, barrier=(40.4315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41377,0.0372982,-2.10076e-05,5.4868e-09,-5.753e-13,20052,16.5196], Tmin=(100,'K'), Tmax=(2046.51,'K')), NASAPolynomial(coeffs=[10.5294,0.0214357,-9.38114e-06,1.69941e-09,-1.12635e-13,16730.2,-28.4458], Tmin=(2046.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]O[C]=O(587)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (31.9823,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,459.134],'cm^-1')),
        HinderedRotor(inertia=(0.110879,'amu*angstrom^2'), symmetry=1, barrier=(16.5839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110914,'amu*angstrom^2'), symmetry=1, barrier=(16.5835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66194,0.0240947,-1.5247e-05,-2.54387e-09,3.90632e-12,3899.55,13.747], Tmin=(100,'K'), Tmax=(979.978,'K')), NASAPolynomial(coeffs=[10.9901,0.0022878,-5.21442e-07,1.28003e-10,-1.24758e-14,1682.1,-29.2488], Tmin=(979.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CsJOC(O)H) + radical((O)CJOCH3)"""),
)

species(
    label = '[CH]C=C(1759)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 C u2 p0 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,227.33,229.304,232.839],'cm^-1')),
        HinderedRotor(inertia=(1.31513,'amu*angstrom^2'), symmetry=1, barrier=(50.4662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31911,0.00817973,3.34731e-05,-4.36187e-08,1.5821e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.759,'K')), NASAPolynomial(coeffs=[5.3676,0.0170742,-6.35103e-06,1.16618e-09,-8.27612e-14,44095,-3.44632], Tmin=(983.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1COC1=O(1882)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-223.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24431,0.0237125,5.22147e-05,-8.56621e-08,3.55018e-11,-26760.4,21.2712], Tmin=(100,'K'), Tmax=(923.422,'K')), NASAPolynomial(coeffs=[11.678,0.0210669,-5.56978e-06,8.75807e-10,-6.10638e-14,-30132.1,-32.3125], Tmin=(923.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Beta-Propiolactone)"""),
)

species(
    label = 'C=CC=COC=O(1875)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,D} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-242.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966509,0.0519695,-9.17413e-06,-3.36432e-08,1.8612e-11,-29034.8,22.6498], Tmin=(100,'K'), Tmax=(966.768,'K')), NASAPolynomial(coeffs=[18.7428,0.0136845,-4.48756e-06,8.55519e-10,-6.59835e-14,-34119.8,-71.0328], Tmin=(966.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH)"""),
)

species(
    label = 'C=C=CCOC=O(1885)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,D} {10,S}
5  C u0 p0 c0 {1,S} {2,D} {11,S}
6  C u0 p0 c0 {7,D} {12,S} {13,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-157.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50631,0.0483603,-2.82417e-05,7.19742e-09,-6.22801e-13,-18845.3,24.4704], Tmin=(100,'K'), Tmax=(1626.85,'K')), NASAPolynomial(coeffs=[14.2066,0.0227743,-9.85163e-06,1.7926e-09,-1.19747e-13,-23724,-45.2756], Tmin=(1626.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O=[C]OCC1[CH]C1(2009)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {4,S} {13,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (121.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65856,0.0405207,7.38624e-07,-2.88843e-08,1.34225e-11,14763.7,25.6764], Tmin=(100,'K'), Tmax=(1012.65,'K')), NASAPolynomial(coeffs=[12.6662,0.0219732,-8.71954e-06,1.65583e-09,-1.19522e-13,11255.9,-33.8796], Tmin=(1012.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdOsH) + ring(Cyclopropane) + radical(cyclopropane) + radical((O)CJOCC)"""),
)

species(
    label = 'O=C1C[CH][CH]CO1(2010)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {6,S} {13,S}
6  C u1 p0 c0 {4,S} {5,S} {12,S}
7  C u0 p0 c0 {1,S} {2,D} {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-34.1199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5111,0.0128463,8.39242e-05,-1.09745e-07,3.94338e-11,-4032.68,21.5509], Tmin=(100,'K'), Tmax=(1019.8,'K')), NASAPolynomial(coeffs=[11.9502,0.0271117,-1.24983e-05,2.60633e-09,-1.98633e-13,-8624.83,-37.2481], Tmin=(1019.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCJCC=O) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]COC1=O(2011)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {4,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u1 p0 c0 {3,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1.23237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19411,0.0235967,5.62203e-05,-9.27268e-08,3.88191e-11,-68.2467,24.1432], Tmin=(100,'K'), Tmax=(918.103,'K')), NASAPolynomial(coeffs=[12.7781,0.0190886,-4.38737e-06,6.40131e-10,-4.48337e-14,-3765.12,-35.5633], Tmin=(918.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.23237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(butyrolactone) + radical(CCJCO) + radical(CJC(C)C=O)"""),
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
    label = '[O][C]=O(803)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186753,1.03202e-05,-1.52373e-08,5.80537e-12,3804.5,8.40407], Tmin=(100,'K'), Tmax=(1021.26,'K')), NASAPolynomial(coeffs=[6.36178,0.000422737,-4.06602e-07,1.525e-10,-1.51974e-14,2816.75,-6.43921], Tmin=(1021.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = 'C=CC=CO[C]=O(2012)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,D} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-45.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27608,'amu*angstrom^2'), symmetry=1, barrier=(29.3396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2773,'amu*angstrom^2'), symmetry=1, barrier=(29.3677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27567,'amu*angstrom^2'), symmetry=1, barrier=(29.3302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0919,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763568,0.0596597,-3.70324e-05,-5.66994e-09,9.38137e-12,-5403.01,22.0055], Tmin=(100,'K'), Tmax=(956.98,'K')), NASAPolynomial(coeffs=[19.0167,0.0106976,-3.12963e-06,5.57422e-10,-4.21916e-14,-10148.2,-71.7918], Tmin=(956.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = 'C=C=CCO[C]=O(2013)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {6,D} {11,S} {12,S}
6  C u0 p0 c0 {4,D} {5,D}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (40.2015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,180,369.765,1424.43],'cm^-1')),
        HinderedRotor(inertia=(0.777495,'amu*angstrom^2'), symmetry=1, barrier=(17.8761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777511,'amu*angstrom^2'), symmetry=1, barrier=(17.8765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.77752,'amu*angstrom^2'), symmetry=1, barrier=(17.8767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0919,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47742,0.0528709,-4.42342e-05,1.90988e-08,-3.35757e-12,4927.97,24.1781], Tmin=(100,'K'), Tmax=(1339.98,'K')), NASAPolynomial(coeffs=[11.9103,0.0217274,-9.37173e-06,1.75407e-09,-1.21575e-13,2131.99,-29.2079], Tmin=(1339.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.2015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical((O)CJOCC)"""),
)

species(
    label = 'C=CC[CH]O[C]=O(1997)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u1 p0 c0 {1,S} {3,S} {11,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (64.9641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,357.67,357.673,357.676,357.678],'cm^-1')),
        HinderedRotor(inertia=(0.00131774,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165212,'amu*angstrom^2'), symmetry=1, barrier=(14.998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165207,'amu*angstrom^2'), symmetry=1, barrier=(14.998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165213,'amu*angstrom^2'), symmetry=1, barrier=(14.998,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600864,0.0667805,-6.58852e-05,3.25332e-08,-6.2689e-12,7942.49,27.8775], Tmin=(100,'K'), Tmax=(1270.9,'K')), NASAPolynomial(coeffs=[16.7776,0.0158661,-5.79256e-06,1.01073e-09,-6.80791e-14,3830.68,-54.044], Tmin=(1270.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.9641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(CCsJOC(O)H) + radical((O)CJOCC)"""),
)

species(
    label = 'C=[C]CCO[C]=O(2014)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {5,D}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (108.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180,409.624,409.625,409.625],'cm^-1')),
        HinderedRotor(inertia=(0.0179934,'amu*angstrom^2'), symmetry=1, barrier=(2.14245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0179933,'amu*angstrom^2'), symmetry=1, barrier=(2.14245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709361,'amu*angstrom^2'), symmetry=1, barrier=(16.3096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136976,'amu*angstrom^2'), symmetry=1, barrier=(16.3096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37847,0.0575334,-5.03768e-05,2.3507e-08,-4.53726e-12,13099.1,26.0083], Tmin=(100,'K'), Tmax=(1216.59,'K')), NASAPolynomial(coeffs=[10.9869,0.0259417,-1.14253e-05,2.16224e-09,-1.51038e-13,10761.2,-22.2307], Tmin=(1216.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_S) + radical((O)CJOCC)"""),
)

species(
    label = '[CH2]C=C[CH]OC=O(2015)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,D} {5,S} {9,S}
4  C u0 p0 c0 {3,D} {6,S} {8,S}
5  C u1 p0 c0 {1,S} {3,S} {10,S}
6  C u1 p0 c0 {4,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,D} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-71.6463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70799,0.0436343,-1.52772e-05,-3.87617e-09,2.46991e-12,-8529.11,25.899], Tmin=(100,'K'), Tmax=(1311.55,'K')), NASAPolynomial(coeffs=[10.8803,0.0286079,-1.28998e-05,2.44235e-09,-1.69244e-13,-12048.7,-25.0851], Tmin=(1311.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.6463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(C=CCJ(O)C) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CCCO[C]=O(2016)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {12,S}
6  C u1 p0 c0 {5,D} {13,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (117.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1855,455,950,269.827,269.83,269.833],'cm^-1')),
        HinderedRotor(inertia=(0.00231531,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369657,'amu*angstrom^2'), symmetry=1, barrier=(19.0993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369654,'amu*angstrom^2'), symmetry=1, barrier=(19.0993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00231531,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12025,0.0589893,-5.09308e-05,2.2597e-08,-4.03848e-12,14225,26.8658], Tmin=(100,'K'), Tmax=(1332.86,'K')), NASAPolynomial(coeffs=[13.4995,0.0218386,-9.12151e-06,1.685e-09,-1.16112e-13,10925,-36.4138], Tmin=(1332.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_P) + radical((O)CJOCC)"""),
)

species(
    label = 'C=[C][CH]COC=O(2017)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u1 p0 c0 {3,S} {7,S} {10,S}
5  C u0 p0 c0 {1,S} {2,D} {13,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (27.3623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,395.563,395.564,395.564,395.568],'cm^-1')),
        HinderedRotor(inertia=(0.00107734,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398179,'amu*angstrom^2'), symmetry=1, barrier=(44.2121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398177,'amu*angstrom^2'), symmetry=1, barrier=(44.212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398175,'amu*angstrom^2'), symmetry=1, barrier=(44.2121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31004,0.0558969,-3.96819e-05,1.37129e-08,-1.92033e-12,3389.98,23.822], Tmin=(100,'K'), Tmax=(1633.22,'K')), NASAPolynomial(coeffs=[13.9402,0.0249636,-1.12716e-05,2.11592e-09,-1.45153e-13,-735.543,-43.3067], Tmin=(1633.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.3623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH]C=CCOC=O(2018)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {12,S}
7  C u2 p0 c0 {5,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (36.5998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48457,0.0485312,-1.87834e-05,-1.49907e-09,1.74457e-12,4497.95,26.5957], Tmin=(100,'K'), Tmax=(1362.18,'K')), NASAPolynomial(coeffs=[10.9245,0.0326021,-1.42266e-05,2.62532e-09,-1.78634e-13,832.261,-25.8796], Tmin=(1362.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.5998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(AllylJ2_triplet)"""),
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
    label = 'C=C[CH]O[C]=O(1264)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u1 p0 c0 {1,S} {3,S} {8,S}
5  C u0 p0 c0 {3,D} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {2,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (10.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,379.487,379.488,2036.76],'cm^-1')),
        HinderedRotor(inertia=(0.35913,'amu*angstrom^2'), symmetry=1, barrier=(36.701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0067017,'amu*angstrom^2'), symmetry=1, barrier=(19.7282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359133,'amu*angstrom^2'), symmetry=1, barrier=(36.7011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1581,0.0367455,-2.39004e-05,7.48573e-09,-9.41401e-13,1339.58,21.8223], Tmin=(100,'K'), Tmax=(1816.37,'K')), NASAPolynomial(coeffs=[11.4939,0.016186,-6.92163e-06,1.25388e-09,-8.36541e-14,-2051.83,-28.7895], Tmin=(1816.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=CCJ(O)C) + radical((O)CJOCC)"""),
)

species(
    label = '[CH2]C(C=C)O[C]=O(1881)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (65.3836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,350.995,351.706,351.762],'cm^-1')),
        HinderedRotor(inertia=(0.00136822,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179211,'amu*angstrom^2'), symmetry=1, barrier=(15.6669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179291,'amu*angstrom^2'), symmetry=1, barrier=(15.6658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179134,'amu*angstrom^2'), symmetry=1, barrier=(15.6649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3602.61,'J/mol'), sigma=(6.00852,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.72 K, Pc=37.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885483,0.0652967,-6.43073e-05,3.24103e-08,-6.50333e-12,7978.66,27.3716], Tmin=(100,'K'), Tmax=(1204.67,'K')), NASAPolynomial(coeffs=[14.3456,0.0206037,-8.65771e-06,1.6138e-09,-1.1227e-13,4735.65,-40.0721], Tmin=(1204.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.3836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(CJC(C)OC) + radical((O)CJOCC2)"""),
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
    label = '[CH]=CCO[C]=O(1485)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u1 p0 c0 {4,D} {10,S}
6  C u1 p0 c0 {1,S} {2,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (146.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1855,455,950,240.037,240.047],'cm^-1')),
        HinderedRotor(inertia=(0.0997365,'amu*angstrom^2'), symmetry=1, barrier=(4.07808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.448586,'amu*angstrom^2'), symmetry=1, barrier=(18.3418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0171813,'amu*angstrom^2'), symmetry=1, barrier=(18.3419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89904,0.0424688,-3.57286e-05,1.53434e-08,-2.65022e-12,17724.9,22.2101], Tmin=(100,'K'), Tmax=(1377.11,'K')), NASAPolynomial(coeffs=[11.134,0.0156443,-6.51007e-06,1.19838e-09,-8.23078e-14,15181.4,-25.2981], Tmin=(1377.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_P) + radical((O)CJOCC)"""),
)

species(
    label = 'O=C1CC=CCO1(1884)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {12,S}
6  C u0 p0 c0 {4,S} {5,D} {13,S}
7  C u0 p0 c0 {1,S} {2,D} {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-310.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04166,-0.00156053,0.000123926,-1.44079e-07,4.90507e-11,-37290,14.7904], Tmin=(100,'K'), Tmax=(1038.26,'K')), NASAPolynomial(coeffs=[9.70799,0.034383,-1.70347e-05,3.59942e-09,-2.7365e-13,-41995.9,-33.6172], Tmin=(1038.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-310.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclohexane)"""),
)

species(
    label = '[CH2][CH]C1COC1=O(2019)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u1 p0 c0 {5,S} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (55.3625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22231,0.0234223,5.47372e-05,-9.04471e-08,3.78808e-11,6737.12,25.8272], Tmin=(100,'K'), Tmax=(917.457,'K')), NASAPolynomial(coeffs=[12.4838,0.0191567,-4.46055e-06,6.52217e-10,-4.53927e-14,3150.85,-32.0777], Tmin=(917.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.3625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = 'C[C]=CCO[C]=O(2020)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u1 p0 c0 {4,S} {5,D}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (101.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6739,0.0545659,-4.82734e-05,2.45687e-08,-5.43834e-12,12281.2,24.9014], Tmin=(100,'K'), Tmax=(1040.72,'K')), NASAPolynomial(coeffs=[7.82764,0.030914,-1.41837e-05,2.73148e-09,-1.92638e-13,11000.3,-5.0325], Tmin=(1040.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Cds_S) + radical((O)CJOCC)"""),
)

species(
    label = 'CC=[C]CO[C]=O(2021)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,D} {13,S}
6  C u1 p0 c0 {3,S} {5,D}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (101.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,458.745,458.746],'cm^-1')),
        HinderedRotor(inertia=(0.0142331,'amu*angstrom^2'), symmetry=1, barrier=(2.12553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0687561,'amu*angstrom^2'), symmetry=1, barrier=(10.2679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.446586,'amu*angstrom^2'), symmetry=1, barrier=(10.2679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.068756,'amu*angstrom^2'), symmetry=1, barrier=(10.2679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6739,0.0545659,-4.82734e-05,2.45687e-08,-5.43834e-12,12281.2,24.9014], Tmin=(100,'K'), Tmax=(1040.72,'K')), NASAPolynomial(coeffs=[7.82764,0.030914,-1.41837e-05,2.73148e-09,-1.92638e-13,11000.3,-5.0325], Tmin=(1040.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Cds_S) + radical((O)CJOCC)"""),
)

species(
    label = '[CH2]C=[C]COC=O(2022)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {7,D} {10,S}
5  C u1 p0 c0 {4,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {2,D} {13,S}
7  C u1 p0 c0 {3,S} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (55.2563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1685,370,377.196,377.197,1795.75],'cm^-1')),
        HinderedRotor(inertia=(0.373268,'amu*angstrom^2'), symmetry=1, barrier=(37.6864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0921416,'amu*angstrom^2'), symmetry=1, barrier=(9.30284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373265,'amu*angstrom^2'), symmetry=1, barrier=(37.6864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164689,'amu*angstrom^2'), symmetry=1, barrier=(37.6864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59529,0.0491101,-3.01881e-05,8.71556e-09,-1.00255e-12,6735.18,25.8367], Tmin=(100,'K'), Tmax=(1966.66,'K')), NASAPolynomial(coeffs=[14.8856,0.0220791,-9.57121e-06,1.72681e-09,-1.14154e-13,1507.65,-47.2701], Tmin=(1966.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.2563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[CH]O[C]=O(2023)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u1 p0 c0 {1,S} {5,S} {13,S}
7  C u1 p0 c0 {1,S} {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-25.4636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88861,0.0481848,-3.12072e-05,1.00528e-08,-1.35634e-12,-2988.78,24.5723], Tmin=(100,'K'), Tmax=(1607.75,'K')), NASAPolynomial(coeffs=[9.85959,0.0283534,-1.27049e-05,2.38061e-09,-1.63335e-13,-5551.85,-17.6679], Tmin=(1607.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.4636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + radical(C=CCJ(O)C) + radical((O)CJOCC)"""),
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
    E0 = (-78.3901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (539.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (343.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-70.1057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-0.142978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-53.4168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (147.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-18.6766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-2.16814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-4.8684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-63.7941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (73.6165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (101.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (186.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (240.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (106.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (200.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (86.3918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (196.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (113.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-3.85511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (515.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (192.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (462.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-70.8589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-3.74367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (129.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (197.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (90.8414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (39.1803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['CO2(15)', 'C=CC=C(61)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=O(152)', 'C=C[CH]C[O](1223)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]O[C]=O(587)', '[CH]C=C(1759)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['C=CC1COC1=O(1882)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['C=CC=COC=O(1875)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['C=C=CCOC=O(1885)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['O=[C]OCC1[CH]C1(2009)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['O=C1C[CH][CH]CO1(2010)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.25216e+09,'s^-1'), n=0.38532, Ea=(59.7134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['[CH2]C1[CH]COC1=O(2011)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.70833e+09,'s^-1'), n=0.576561, Ea=(76.2219,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CO(14)', 'C=C[CH]C[O](1223)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(34.1,'m^3/(mol*s)'), n=8.73864e-09, Ea=(13.1922,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO2(15)', '[CH2]C=C[CH2](74)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.64086e-07,'m^3/(mol*s)'), n=3.71185, Ea=(130.222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(803)', 'C=CC=C(61)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.64086e-07,'m^3/(mol*s)'), n=3.71185, Ea=(11.2182,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(6)', 'C=CC=CO[C]=O(2012)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(21820,'m^3/(mol*s)'), n=0.859, Ea=(0.819436,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(6)', 'C=C=CCO[C]=O(2013)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.35489,'m^3/(mol*s)'), n=2.14763, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.22132624170300447, var=0.24654008949797426, Tref=1000.0, N=22, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_N-4R!H->O_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C]=O(803)', '[CH2]C=C[CH2](74)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.95326e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=CC[CH]O[C]=O(1997)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.46162e+08,'s^-1'), n=1.28739, Ea=(107.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]CCO[C]=O(2014)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['[CH2]C=C[CH]OC=O(2015)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.8344e+08,'s^-1'), n=1.32036, Ea=(164.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;CO_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CCCO[C]=O(2016)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C][CH]COC=O(2017)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C=CCOC=O(2018)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.546e+09,'s^-1'), n=0.732, Ea=(25.1375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(24)', 'C=C[CH]O[C]=O(1264)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(151.149,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)O[C]=O(1881)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(192.885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(17)', '[CH]=CCO[C]=O(1485)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['O=C1CC=CCO1(1884)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['[CH2][CH]C1COC1=O(2019)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(74.6464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['C[C]=CCO[C]=O(2020)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC=[C]CO[C]=O(2021)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C=[C]COC=O(2022)'],
    products = ['C=C[CH]CO[C]=O(1880)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[CH]CO[C]=O(1880)'],
    products = ['CC=C[CH]O[C]=O(2023)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(256000,'s^-1'), n=2, Ea=(117.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_2H;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1116',
    isomers = [
        'C=C[CH]CO[C]=O(1880)',
    ],
    reactants = [
        ('CO2(15)', 'C=CC=C(61)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1116',
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

