species(
    label = '[CH]=C(C)C([O])=O(1174)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {5,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {1,S} {2,D} {4,S}
6  C u1 p0 c0 {4,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (175.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,180,180,180,1332.24,3627.54],'cm^-1')),
        HinderedRotor(inertia=(0.149379,'amu*angstrom^2'), symmetry=1, barrier=(3.43451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00272671,'amu*angstrom^2'), symmetry=1, barrier=(3.43464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.73068,0.0245579,-9.08827e-06,2.88061e-10,1.65859e-13,21053.1,11.0378], Tmin=(100,'K'), Tmax=(2729.23,'K')), NASAPolynomial(coeffs=[41.5047,-0.0130867,2.66928e-06,-4.02093e-10,2.92172e-14,-5072.64,-214.386], Tmin=(2729.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
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
    label = '[CH]=CC([O])=O(851)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 C u1 p0 c0 {3,D} {7,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (185.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000717289,'amu*angstrom^2'), symmetry=1, barrier=(8.14414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3904.09,'J/mol'), sigma=(5.73968,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=609.81 K, Pc=46.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36226,0.0216394,-5.09964e-05,6.30528e-08,-2.69918e-11,22291.4,13.025], Tmin=(100,'K'), Tmax=(834.91,'K')), NASAPolynomial(coeffs=[-0.323966,0.0192643,-1.07333e-05,2.16075e-09,-1.52238e-13,23605.2,34.3258], Tmin=(834.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=C([O])OC(1268)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {6,D}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (131.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29287,0.0556112,-5.89564e-05,2.86934e-08,-4.44901e-12,15907.5,19.2537], Tmin=(100,'K'), Tmax=(927.132,'K')), NASAPolynomial(coeffs=[13.9485,0.0111422,-3.40317e-06,5.34464e-10,-3.44421e-14,13125.4,-43.1932], Tmin=(927.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = 'O(7)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(C)=C=O(1269)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3 C u0 p0 c0 {2,S} {4,S} {5,D}
4 C u2 p0 c0 {3,S} {9,S}
5 C u0 p0 c0 {1,D} {3,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (267.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2120,512.5,787.5,402.449,402.451,402.455],'cm^-1')),
        HinderedRotor(inertia=(0.453263,'amu*angstrom^2'), symmetry=1, barrier=(52.0957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453257,'amu*angstrom^2'), symmetry=1, barrier=(52.0958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.0738,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60797,0.0336788,-1.93865e-05,5.97209e-09,-8.55154e-13,32210,15.5845], Tmin=(100,'K'), Tmax=(1346.14,'K')), NASAPolynomial(coeffs=[4.92928,0.0267811,-1.17004e-05,2.16562e-09,-1.4823e-13,31585.1,3.69559], Tmin=(1346.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + missing(Cdd-CdO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CC1=COC1=O(1261)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {1,S} {4,D} {10,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-122.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06824,0.0357555,-1.71571e-05,-2.89681e-09,3.31668e-12,-14705.8,16.5914], Tmin=(100,'K'), Tmax=(1086.71,'K')), NASAPolynomial(coeffs=[10.9809,0.0155239,-6.58775e-06,1.26699e-09,-9.07545e-14,-17385.4,-30.5645], Tmin=(1086.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C(C)[C]1OO1(1270)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  C u1 p0 c0 {4,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (453.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76266,0.0410407,-1.74392e-05,-1.47439e-08,1.12913e-11,54639.3,20.4562], Tmin=(100,'K'), Tmax=(912.307,'K')), NASAPolynomial(coeffs=[13.8714,0.00980921,-2.02929e-06,2.58705e-10,-1.68655e-14,51520.2,-41.8356], Tmin=(912.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = 'C[C]1[CH]OC1=O(1271)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  C u1 p0 c0 {1,S} {4,S} {10,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (20.4956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61537,0.0149605,5.8246e-05,-9.31012e-08,3.9524e-11,2529.43,18.2519], Tmin=(100,'K'), Tmax=(903.141,'K')), NASAPolynomial(coeffs=[13.0052,0.00955664,-2.30786e-07,-1.44931e-10,9.95455e-15,-1003.56,-39.9841], Tmin=(903.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.4956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCJ(C)CO) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC1=CC1([O])[O](1272)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (311.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02169,0.0499857,-7.45173e-05,7.25491e-08,-2.83915e-11,37517.7,17.2778], Tmin=(100,'K'), Tmax=(799.86,'K')), NASAPolynomial(coeffs=[2.9422,0.0321174,-1.61323e-05,3.15262e-09,-2.20889e-13,37794.8,15.695], Tmin=(799.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'C#CC([O])=O(941)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u0 p0 c0 {4,T} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (137.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180,180,1347.07,1348.7,1352.25],'cm^-1')),
        HinderedRotor(inertia=(0.242715,'amu*angstrom^2'), symmetry=1, barrier=(5.58049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0387,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85456,0.034387,-7.33086e-05,7.96015e-08,-3.07521e-11,16569.1,14.0024], Tmin=(100,'K'), Tmax=(878.45,'K')), NASAPolynomial(coeffs=[0.760991,0.0207821,-1.05685e-05,2.00329e-09,-1.35032e-13,17829.7,28.9127], Tmin=(878.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ)"""),
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
    label = '[CH2]C(=C)C([O])=O(1273)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u1 p0 c0 {3,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (84.8686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87128,0.0287993,-1.40179e-05,2.55385e-09,-1.57424e-13,10183.9,12.5909], Tmin=(100,'K'), Tmax=(3036.44,'K')), NASAPolynomial(coeffs=[40.7824,-0.0129792,3.23888e-06,-4.92499e-10,3.2261e-14,-15387.6,-211.677], Tmin=(3036.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.8686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH]C(=C)C(=O)O(1274)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u0 p0 c0 {1,S} {2,D} {3,S}
5  C u0 p0 c0 {3,D} {7,S} {8,S}
6  C u2 p0 c0 {3,S} {9,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
"""),
    E0 = (73.4965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92159,0.0505504,-6.34918e-05,5.63118e-08,-2.13755e-11,8909.74,21.8948], Tmin=(100,'K'), Tmax=(779.085,'K')), NASAPolynomial(coeffs=[3.68848,0.0337416,-1.62326e-05,3.12472e-09,-2.17977e-13,8869.24,15.3186], Tmin=(779.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.4965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (18.3624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (450.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (387.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (353.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (26.6467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (296.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (156.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (154.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (148.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (103.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (76.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (365.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (147.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (18.3624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['CO2(15)', 'C#CC(211)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(24)', '[CH]=CC([O])=O(851)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.11579e+08,'m^3/(mol*s)'), n=-0.607843, Ea=(3.45865,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['[CH]=C=C([O])OC(1268)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6862.26,'s^-1'), n=2.77015, Ea=(368.922,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(7)', '[CH]C(C)=C=O(1269)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['CC1=COC1=O(1261)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['[CH]=C(C)[C]1OO1(1270)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(277.95,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_De;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 277.2 to 278.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['C[C]1[CH]OC1=O(1271)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.906e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['CC1=CC1([O])[O](1272)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(135.844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(18)', 'C#CC([O])=O(941)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00276756,'m^3/(mol*s)'), n=2.82919, Ea=(32.308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.07604778385317046, var=1.4951309799213957, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=O(803)', 'C#CC(211)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.19389,'m^3/(mol*s)'), n=1.76036, Ea=(56.754,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.10151511122778119, var=2.139043046412222, Tref=1000.0, N=264, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO2(15)', '[CH]=[C]C(212)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(145.878,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(803)', '[CH]=[C]C(212)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C)C([O])=O(1174)'],
    products = ['[CH2]C(=C)C([O])=O(1273)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(=C)C(=O)O(1274)'],
    products = ['[CH]=C(C)C([O])=O(1174)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(102.117,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 101.5 to 102.1 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #709',
    isomers = [
        '[CH]=C(C)C([O])=O(1174)',
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
    label = 'PDepNetwork #709',
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

