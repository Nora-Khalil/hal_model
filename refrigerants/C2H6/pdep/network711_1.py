species(
    label = '[CH]=C(C)O[C]=O(1176)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u1 p0 c0 {4,D} {10,S}
6  C u1 p0 c0 {1,S} {2,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (107.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.867935,'amu*angstrom^2'), symmetry=1, barrier=(19.9555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867921,'amu*angstrom^2'), symmetry=1, barrier=(19.9552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867846,'amu*angstrom^2'), symmetry=1, barrier=(19.9535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68062,0.0535051,-6.07889e-05,3.66561e-08,-8.98498e-12,13017.2,18.4957], Tmin=(100,'K'), Tmax=(982.64,'K')), NASAPolynomial(coeffs=[9.93474,0.0199058,-9.50028e-06,1.86007e-09,-1.32417e-13,11395,-21.1812], Tmin=(982.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_P) + radical((O)CJOC)"""),
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
    label = 'C#CC(211)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (172.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,743.639],'cm^-1')),
        HinderedRotor(inertia=(0.630823,'amu*angstrom^2'), symmetry=1, barrier=(14.5039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30524,0.0109264,1.31992e-05,-2.25168e-08,8.87439e-12,20738.2,7.1916], Tmin=(100,'K'), Tmax=(969.981,'K')), NASAPolynomial(coeffs=[5.80052,0.0113536,-4.03494e-06,7.1907e-10,-5.02191e-14,19750,-7.36949], Tmin=(969.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=CO[C]=O(852)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,D} {6,S}
4 C u1 p0 c0 {3,D} {7,S}
5 C u1 p0 c0 {1,S} {2,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (149.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.44717,'amu*angstrom^2'), symmetry=1, barrier=(33.2732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44585,'amu*angstrom^2'), symmetry=1, barrier=(33.243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3289.41,'J/mol'), sigma=(5.34165,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.80 K, Pc=48.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72784,0.0416654,-4.5885e-05,2.33989e-08,-4.46912e-12,18051,16.3166], Tmin=(100,'K'), Tmax=(1415.99,'K')), NASAPolynomial(coeffs=[14.2323,0.00248459,-2.93512e-07,1.00895e-11,-5.60316e-17,14896.4,-46.9941], Tmin=(1415.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_P) + radical((O)CJOC)"""),
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
    label = '[CH]=C(C)[O](245)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u1 p0 c0 {3,D} {8,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (198.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.320495,'amu*angstrom^2'), symmetry=1, barrier=(7.36881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0632,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64262,0.0274443,-2.15887e-05,9.33787e-09,-1.65053e-12,23935.7,14.2603], Tmin=(100,'K'), Tmax=(1341.74,'K')), NASAPolynomial(coeffs=[7.83006,0.0119795,-4.29974e-06,7.47548e-10,-4.99364e-14,22543.6,-12.291], Tmin=(1341.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'CC1=CC(=O)O1(1237)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {6,S} {10,S}
6  C u0 p0 c0 {1,S} {2,D} {5,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-155.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76059,0.0226311,-6.7106e-06,-4.91374e-09,2.57234e-12,-18604.6,12.9591], Tmin=(100,'K'), Tmax=(1192.51,'K')), NASAPolynomial(coeffs=[8.5947,0.0124573,-5.73134e-06,1.14551e-09,-8.29614e-14,-20664.1,-19.0153], Tmin=(1192.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cyclobutene)"""),
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
    label = '[CH]=[C]C(212)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (490.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0859242,'amu*angstrom^2'), symmetry=1, barrier=(1.97557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22084,0.0161542,-7.04072e-06,1.20991e-09,-3.01638e-14,59071.4,10.7688], Tmin=(100,'K'), Tmax=(1807.72,'K')), NASAPolynomial(coeffs=[6.96329,0.0100049,-3.70699e-06,6.32786e-10,-4.05637e-14,57370,-10.4656], Tmin=(1807.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
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
    label = 'CH3(18)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C#CO[C]=O(253)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {5,T}
4 C u1 p0 c0 {1,S} {2,D}
5 C u0 p0 c0 {3,T} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (131.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,1855,455,950,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.83092,'amu*angstrom^2'), symmetry=1, barrier=(42.0964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8282,'amu*angstrom^2'), symmetry=1, barrier=(42.0339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0387,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3118,0.0326921,-3.23756e-05,1.50433e-08,-2.70077e-12,15845.5,12.7313], Tmin=(100,'K'), Tmax=(1362,'K')), NASAPolynomial(coeffs=[11.4718,0.00579011,-2.74751e-06,5.40775e-10,-3.87436e-14,13350.3,-34.2902], Tmin=(1362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-OdOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C(=C)O[C]=O(1238)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {7,S} {8,S}
6  C u1 p0 c0 {1,S} {2,D}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (19.3742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52066,0.052959,-5.62801e-05,3.05731e-08,-6.58404e-12,2420.83,18.6555], Tmin=(100,'K'), Tmax=(1128.87,'K')), NASAPolynomial(coeffs=[12.0737,0.0155656,-6.59288e-06,1.22967e-09,-8.55842e-14,38.2473,-33.5358], Tmin=(1128.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.3742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=C(O)CJ) + radical((O)CJOC)"""),
)

species(
    label = '[CH]C(=C)OC=O(1239)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {3,D} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {2,D} {9,S}
6  C u2 p0 c0 {3,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (34.6935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49861,0.0473816,-3.2784e-05,1.09209e-08,-1.45135e-12,4269.13,21.2425], Tmin=(100,'K'), Tmax=(1754.46,'K')), NASAPolynomial(coeffs=[14.3843,0.0180037,-7.66713e-06,1.37692e-09,-9.14031e-14,-252.382,-48.1674], Tmin=(1754.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(AllylJ2_triplet)"""),
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
    E0 = (-22.4563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (438.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (507.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-14.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1.21938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (23.3486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (94.8265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (161.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (392.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (106.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (79.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C)O[C]=O(1176)'],
    products = ['CO2(15)', 'C#CC(211)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH]=CO[C]=O(852)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.11579e+08,'m^3/(mol*s)'), n=-0.607843, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(152)', '[CH]=C(C)[O](245)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C)O[C]=O(1176)'],
    products = ['CC1=CC(=O)O1(1237)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CO(14)', '[CH]=C(C)[O](245)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(34.1,'m^3/(mol*s)'), n=8.73864e-09, Ea=(51.3805,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CO2(15)', '[CH]=[C]C(212)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0308,'m^3/(mol*s)'), n=2.41, Ea=(65.5893,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_N-Sp-5R!H-2C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R_N-Sp-5R!H-2C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C]=O(803)', 'C#CC(211)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(21.1142,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(18)', 'C#CO[C]=O(253)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(680.244,'m^3/(mol*s)'), n=1.05459, Ea=(24.2132,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3424634175026916, var=2.6292686492138895, Tref=1000.0, N=81, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=O(803)', '[CH]=[C]C(212)'],
    products = ['[CH]=C(C)O[C]=O(1176)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C)O[C]=O(1176)'],
    products = ['[CH2]C(=C)O[C]=O(1238)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C)O[C]=O(1176)'],
    products = ['[CH]C(=C)OC=O(1239)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.52029e+07,'s^-1'), n=1.31708, Ea=(101.589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #711',
    isomers = [
        '[CH]=C(C)O[C]=O(1176)',
    ],
    reactants = [
        ('CO2(15)', 'C#CC(211)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #711',
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

