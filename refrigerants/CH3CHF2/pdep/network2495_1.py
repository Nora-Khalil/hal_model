species(
    label = '[CH]=C(C)C([O])=O(6326)',
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
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,180,180,180,1332.29,3627.52],'cm^-1')),
        HinderedRotor(inertia=(0.149349,'amu*angstrom^2'), symmetry=1, barrier=(3.43383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00272747,'amu*angstrom^2'), symmetry=1, barrier=(3.43544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.73067,0.0245579,-9.08828e-06,2.88065e-10,1.65858e-13,21053.1,11.0378], Tmin=(100,'K'), Tmax=(2729.23,'K')), NASAPolynomial(coeffs=[41.5052,-0.0130873,2.66953e-06,-4.02137e-10,2.92201e-14,-5073.01,-214.39], Tmin=(2729.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
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
    label = 'C#CC(5603)',
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
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,743.64],'cm^-1')),
        HinderedRotor(inertia=(0.630825,'amu*angstrom^2'), symmetry=1, barrier=(14.5039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30522,0.0109266,1.31984e-05,-2.25157e-08,8.87389e-12,20738.2,7.19167], Tmin=(100,'K'), Tmax=(969.992,'K')), NASAPolynomial(coeffs=[5.80059,0.0113535,-4.03487e-06,7.19054e-10,-5.02177e-14,19749.9,-7.36988], Tmin=(969.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144066,5.45064e-06,-3.57997e-09,7.56176e-13,50400.6,-0.411759], Tmin=(100,'K'), Tmax=(1442.38,'K')), NASAPolynomial(coeffs=[2.62651,0.00394757,-1.49921e-06,2.54533e-10,-1.62951e-14,50691.7,6.78356], Tmin=(1442.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=CC([O])=O(683)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36235,0.0216383,-5.09919e-05,6.30461e-08,-2.69884e-11,22291.4,13.0247], Tmin=(100,'K'), Tmax=(834.93,'K')), NASAPolynomial(coeffs=[-0.324133,0.0192646,-1.07335e-05,2.1608e-09,-1.52241e-13,23605.3,34.3267], Tmin=(834.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=C([O])OC(6630)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29287,0.0556112,-5.89567e-05,2.86937e-08,-4.44917e-12,15907.5,19.2538], Tmin=(100,'K'), Tmax=(927.135,'K')), NASAPolynomial(coeffs=[13.9485,0.0111422,-3.40316e-06,5.3446e-10,-3.44418e-14,13125.4,-43.1933], Tmin=(927.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,29230.2,5.12617], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(C)=C=O(6514)',
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
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2120,512.5,787.5,397.905,403.793,405.824],'cm^-1')),
        HinderedRotor(inertia=(0.466414,'amu*angstrom^2'), symmetry=1, barrier=(52.0313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441442,'amu*angstrom^2'), symmetry=1, barrier=(52.1859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.0738,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60797,0.0336788,-1.93865e-05,5.97209e-09,-8.55156e-13,32210,15.5845], Tmin=(100,'K'), Tmax=(1346.14,'K')), NASAPolynomial(coeffs=[4.92927,0.0267811,-1.17004e-05,2.16562e-09,-1.4823e-13,31585.1,3.69561], Tmin=(1346.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + missing(Cdd-CdO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CC1=COC1=O(6623)',
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
    E0 = (-111.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22949,0.0347802,-2.33854e-05,7.58774e-09,-9.82689e-13,-13330.7,15.8501], Tmin=(100,'K'), Tmax=(1785.39,'K')), NASAPolynomial(coeffs=[11.4667,0.0140851,-5.99839e-06,1.09541e-09,-7.35999e-14,-16629.1,-34.0684], Tmin=(1785.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Oxetene)"""),
)

species(
    label = '[CH]=C(C)[C]1OO1(6631)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76266,0.0410407,-1.74392e-05,-1.47439e-08,1.12913e-11,54639.3,20.4562], Tmin=(100,'K'), Tmax=(912.307,'K')), NASAPolynomial(coeffs=[13.8714,0.00980921,-2.02929e-06,2.58705e-10,-1.68654e-14,51520.2,-41.8356], Tmin=(912.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = 'C[C]1[CH]OC1=O(6632)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61535,0.0149608,5.82449e-05,-9.30996e-08,3.95233e-11,2529.43,18.252], Tmin=(100,'K'), Tmax=(903.144,'K')), NASAPolynomial(coeffs=[13.0052,0.00955653,-2.3072e-07,-1.44947e-10,9.95589e-15,-1003.59,-39.9844], Tmin=(903.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.4956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCJ(C)CO) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC1=CC1([O])[O](6633)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02162,0.0499866,-7.45211e-05,7.2555e-08,-2.83946e-11,37517.7,17.278], Tmin=(100,'K'), Tmax=(799.837,'K')), NASAPolynomial(coeffs=[2.94232,0.0321172,-1.61321e-05,3.15259e-09,-2.20886e-13,37794.7,15.6943], Tmin=(799.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'CH3(19)',
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
    label = 'C#CC([O])=O(732)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u0 p0 c0 {4,T} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (87.8995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,575.085,575.089,575.09,575.092,575.092],'cm^-1')),
        HinderedRotor(inertia=(0.000509711,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0387,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42283,0.016444,1.21374e-05,-6.88388e-08,6.82458e-11,10588.7,10.6547], Tmin=(100,'K'), Tmax=(386.352,'K')), NASAPolynomial(coeffs=[3.92521,0.0164858,-8.38135e-06,1.69264e-09,-1.22499e-13,10510.7,8.20233], Tmin=(386.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.8995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: 2-BTP_G4 + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(293)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.00186761,1.03198e-05,-1.52369e-08,5.8052e-12,3804.5,8.40409], Tmin=(100,'K'), Tmax=(1021.27,'K')), NASAPolynomial(coeffs=[6.36181,0.000422681,-4.06569e-07,1.52493e-10,-1.51968e-14,2816.74,-6.4394], Tmin=(1021.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[CH]=[C]C(5604)',
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
        HinderedRotor(inertia=(0.0859241,'amu*angstrom^2'), symmetry=1, barrier=(1.97556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0638,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22084,0.0161542,-7.04067e-06,1.20986e-09,-3.01512e-14,59071.4,10.7688], Tmin=(100,'K'), Tmax=(1807.69,'K')), NASAPolynomial(coeffs=[6.96312,0.0100051,-3.7071e-06,6.32809e-10,-4.05654e-14,57370.1,-10.4646], Tmin=(1807.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=C)C([O])=O(6634)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u1 p0 c0 {3,S} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (84.8686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0732,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87127,0.0287993,-1.40179e-05,2.55385e-09,-1.57425e-13,10183.9,12.5909], Tmin=(100,'K'), Tmax=(3036.46,'K')), NASAPolynomial(coeffs=[40.7838,-0.0129807,3.23948e-06,-4.92603e-10,3.22676e-14,-15388.7,-211.686], Tmin=(3036.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.8686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH]C(=C)C(=O)O(6635)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {9,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u0 p0 c0 {1,S} {2,D} {3,S}
5  C u0 p0 c0 {3,D} {7,S} {8,S}
6  C u2 p0 c0 {3,S} {10,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {6,S}
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92153,0.0505512,-6.34953e-05,5.63173e-08,-2.13785e-11,8909.74,21.895], Tmin=(100,'K'), Tmax=(779.052,'K')), NASAPolynomial(coeffs=[3.68857,0.0337415,-1.62325e-05,3.12469e-09,-2.17974e-13,8869.21,15.3181], Tmin=(779.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.4965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (20.4502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (452.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (389.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (355.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (28.7345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (298.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (158.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (156.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (118.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (105.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (78.4837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (367.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (149.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (20.4502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['CO2(14)', 'C#CC(5603)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', '[CH]=CC([O])=O(683)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.11579e+08,'m^3/(mol*s)'), n=-0.607843, Ea=(3.45865,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R_N-4CbCdCsCt->Cs_Ext-4CbCdCt-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['[CH]=C=C([O])OC(6630)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6862.26,'s^-1'), n=2.77015, Ea=(368.922,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H-inRing_2R!H->C',), comment="""Estimated from node Root_N-1R!H-inRing_2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', '[CH]C(C)=C=O(6514)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['CC1=COC1=O(6623)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['[CH]=C(C)[C]1OO1(6631)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(277.95,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_De;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 277.2 to 278.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['C[C]1[CH]OC1=O(6632)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.906e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['CC1=CC1([O])[O](6633)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(135.844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(19)', 'C#CC([O])=O(732)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00276756,'m^3/(mol*s)'), n=2.82919, Ea=(49.8541,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.07604778385317046, var=1.4951309799213957, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_N-Sp-2R!H=1R!H_N-4R!H-inRing"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=O(293)', 'C#CC(5603)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.19389,'m^3/(mol*s)'), n=1.76036, Ea=(56.754,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.10151511122778119, var=2.139043046412222, Tref=1000.0, N=264, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO2(14)', '[CH]=[C]C(5604)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(145.878,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(293)', '[CH]=[C]C(5604)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C)C([O])=O(6326)'],
    products = ['[CH2]C(=C)C([O])=O(6634)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(=C)C(=O)O(6635)'],
    products = ['[CH]=C(C)C([O])=O(6326)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(102.117,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 101.5 to 102.1 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #2495',
    isomers = [
        '[CH]=C(C)C([O])=O(6326)',
    ],
    reactants = [
        ('CO2(14)', 'C#CC(5603)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2495',
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

