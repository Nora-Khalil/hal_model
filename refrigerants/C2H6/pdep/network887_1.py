species(
    label = '[CH]C([O])=C[C]=O(1578)',
    structure = adjacencyList("""multiplicity 5
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 C u2 p0 c0 {3,S} {8,S}
6 C u1 p0 c0 {2,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (337.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09344,'amu*angstrom^2'), symmetry=1, barrier=(48.1323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0945,'amu*angstrom^2'), symmetry=1, barrier=(48.1567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46965,0.0598019,-8.50277e-05,6.58156e-08,-2.06104e-11,40641.6,21.5237], Tmin=(100,'K'), Tmax=(778.946,'K')), NASAPolynomial(coeffs=[9.00533,0.0211058,-1.05129e-05,2.04269e-09,-1.43112e-13,39467.6,-12.9492], Tmin=(778.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = 'HCCO(20)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (166.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87608,0.0221205,-3.58869e-05,3.05403e-08,-1.01281e-11,20163.4,13.6968], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91479,0.00371409,-1.30137e-06,2.06473e-10,-1.21477e-14,19359.6,-5.50567], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(166.705,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[CH][C]=C[C]=O(519)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,D} {4,S} {6,S}
3 C u1 p0 c0 {2,D} {5,S}
4 C u1 p0 c0 {1,D} {2,S}
5 C u2 p0 c0 {3,S} {7,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (651.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08439,'amu*angstrom^2'), symmetry=1, barrier=(47.9241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08402,'amu*angstrom^2'), symmetry=1, barrier=(47.9156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.0579,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99077,0.0507469,-8.43727e-05,8.05524e-08,-3.01638e-11,78457.4,18.5592], Tmin=(100,'K'), Tmax=(815.101,'K')), NASAPolynomial(coeffs=[4.49755,0.0252843,-1.32952e-05,2.60947e-09,-1.82251e-13,78486,9.65982], Tmin=(815.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=O(179)',
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
    label = '[CH]C(=[CH])[O](3145)',
    structure = adjacencyList("""multiplicity 5
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u2 p0 c0 {2,S} {5,S}
4 C u1 p0 c0 {2,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (569.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,319.245,319.306,319.311],'cm^-1')),
        HinderedRotor(inertia=(0.706987,'amu*angstrom^2'), symmetry=1, barrier=(51.1423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64262,0.0274443,-2.15887e-05,9.33787e-09,-1.65053e-12,68517.9,15.3589], Tmin=(100,'K'), Tmax=(1341.74,'K')), NASAPolynomial(coeffs=[7.83006,0.0119795,-4.29974e-06,7.47548e-10,-4.99364e-14,67125.8,-11.1924], Tmin=(1341.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1[CH]C(=O)O1(3126)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,S} {6,D}
4 C u1 p0 c0 {3,S} {5,S} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u1 p0 c0 {3,D} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (119.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04491,0.0265901,2.83903e-05,-6.76688e-08,3.09411e-11,14424.6,16.586], Tmin=(100,'K'), Tmax=(933.459,'K')), NASAPolynomial(coeffs=[17.6607,0.00156157,1.29962e-06,-2.49049e-10,9.98701e-15,9684.35,-67.4508], Tmin=(933.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(4-Methylene-2-oxetanone) + radical(C=CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(O)=C=C=O(3146)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {8,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,D}
5 C u2 p0 c0 {3,S} {7,S}
6 C u0 p0 c0 {2,D} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (206.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22253,0.0519799,-5.33926e-05,2.68806e-08,-5.17585e-12,24970.8,20.8371], Tmin=(100,'K'), Tmax=(1380.06,'K')), NASAPolynomial(coeffs=[14.8802,0.00873802,-2.4188e-06,3.37032e-10,-1.96931e-14,21549.3,-48.1915], Tmin=(1380.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cdd-(Cdd-O2d)Cds) + missing(Cdd-CddO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C1OC1[C][O](3115)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {3,S} {4,S}
2 O u1 p2 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u2 p0 c0 {2,S} {3,S}
6 C u1 p0 c0 {4,D} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (755.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52884,0.044468,-3.15476e-05,-4.29424e-09,8.25598e-12,91015.5,18.7732], Tmin=(100,'K'), Tmax=(946.323,'K')), NASAPolynomial(coeffs=[17.8951,-0.000426294,1.12102e-06,-1.91455e-10,8.25848e-15,86830.6,-65.0269], Tmin=(946.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCOJ) + radical(CH2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1([O])C=C1[O](3147)',
    structure = adjacencyList("""multiplicity 5
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (645.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72861,0.0414807,-2.78294e-05,-3.29869e-09,6.73742e-12,77776.5,19.8399], Tmin=(100,'K'), Tmax=(954.129,'K')), NASAPolynomial(coeffs=[15.5123,0.00400655,-8.47185e-07,1.59394e-10,-1.45796e-14,74221.7,-50.856], Tmin=(954.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH](3148)',
    structure = adjacencyList("""multiplicity 4
1 C u3 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (585.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (13.0186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.17631,-0.00339738,5.29661e-06,-3.21805e-09,7.28335e-13,70356.4,-0.992396], Tmin=(100,'K'), Tmax=(1260.71,'K')), NASAPolynomial(coeffs=[3.26552,0.000229838,1.03492e-07,-7.93401e-12,-2.40731e-16,70527.4,3.38021], Tmin=(1260.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CJ3)"""),
)

species(
    label = 'O=[C]C=C=O(470)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u1 p0 c0 {2,D} {3,S}
5 C u0 p0 c0 {1,D} {3,D}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-7.87925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1855,455,950,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.841576,'amu*angstrom^2'), symmetry=1, barrier=(19.3495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0387,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98472,0.0259268,-2.97921e-05,2.05264e-08,-6.41258e-12,-914.289,13.446], Tmin=(100,'K'), Tmax=(741.73,'K')), NASAPolynomial(coeffs=[4.88246,0.0156926,-9.0953e-06,1.92406e-09,-1.42625e-13,-1195.81,4.85751], Tmin=(741.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.87925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + missing(Cdd-CdO2d) + radical(CsCJ=O)"""),
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
    label = '[CH]C([O])=C=C=O(3149)',
    structure = adjacencyList("""multiplicity 4
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,D}
5 C u2 p0 c0 {3,S} {7,S}
6 C u0 p0 c0 {2,D} {4,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (344.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,2120,512.5,787.5,185.915,198.224,199.325,209.139],'cm^-1')),
        HinderedRotor(inertia=(1.81698,'amu*angstrom^2'), symmetry=1, barrier=(50.338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0493,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92604,0.0440777,-4.74083e-05,2.68127e-08,-6.03935e-12,41513.1,19.8106], Tmin=(100,'K'), Tmax=(1081.69,'K')), NASAPolynomial(coeffs=[10.1062,0.013828,-5.46029e-06,9.59172e-10,-6.40807e-14,39743.5,-20.2963], Tmin=(1081.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cdd-(Cdd-O2d)Cds) + missing(Cdd-CddO2d) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([O])=[C][C]=O(3150)',
    structure = adjacencyList("""multiplicity 6
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {3,D} {6,S}
5 C u2 p0 c0 {3,S} {7,S}
6 C u1 p0 c0 {2,D} {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (581.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.03167,'amu*angstrom^2'), symmetry=1, barrier=(46.7121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04948,'amu*angstrom^2'), symmetry=1, barrier=(47.1216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0493,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24229,0.0658411,-0.000109375,9.26917e-08,-3.08331e-11,69990.4,23.0381], Tmin=(100,'K'), Tmax=(803.627,'K')), NASAPolynomial(coeffs=[9.98163,0.0172803,-9.28698e-06,1.82455e-09,-1.27175e-13,68749.2,-16.1967], Tmin=(803.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C([O])=C=C[O](3151)',
    structure = adjacencyList("""multiplicity 5
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {5,D} {6,S}
4 C u0 p0 c0 {2,S} {5,D} {7,S}
5 C u0 p0 c0 {3,D} {4,D}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (373.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08969,'amu*angstrom^2'), symmetry=1, barrier=(48.046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683758,0.0557309,-5.65794e-05,2.72969e-08,-4.87557e-12,45032.2,23.8771], Tmin=(100,'K'), Tmax=(1613.4,'K')), NASAPolynomial(coeffs=[16.5016,0.00578834,-1.7464e-07,-1.3056e-10,1.28789e-14,41324.2,-55.6748], Tmin=(1613.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(O)=[C][C]=O(3152)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {3,S} {8,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {3,D} {6,S}
5 C u2 p0 c0 {3,S} {7,S}
6 C u1 p0 c0 {2,D} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (443.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.91693,'amu*angstrom^2'), symmetry=1, barrier=(44.074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91619,'amu*angstrom^2'), symmetry=1, barrier=(44.0569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91722,'amu*angstrom^2'), symmetry=1, barrier=(44.0807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05643,0.0675265,-9.30035e-05,6.26689e-08,-1.65225e-11,53425.7,22.2144], Tmin=(100,'K'), Tmax=(931.496,'K')), NASAPolynomial(coeffs=[13.4148,0.014458,-7.54714e-06,1.5087e-09,-1.08112e-13,51123.3,-36.5306], Tmin=(931.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH][C]=O(474)',
    structure = adjacencyList("""multiplicity 4
1 O u0 p2 c0 {3,D}
2 C u2 p0 c0 {3,S} {4,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (424.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1287.31,1288.39],'cm^-1')),
        HinderedRotor(inertia=(0.195129,'amu*angstrom^2'), symmetry=1, barrier=(4.48641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62586,0.0087374,-5.11465e-06,1.32025e-09,-1.2688e-13,51044.1,10.36], Tmin=(100,'K'), Tmax=(2333.41,'K')), NASAPolynomial(coeffs=[8.04234,0.00158967,-7.91842e-07,1.62913e-10,-1.121e-14,48867.9,-14.9361], Tmin=(2333.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[C]C([O])=C[C]=O(3153)',
    structure = adjacencyList("""multiplicity 6
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,D} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {6,S}
5 C u1 p0 c0 {2,D} {3,S}
6 C u3 p0 c0 {4,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (635.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.894193,'amu*angstrom^2'), symmetry=1, barrier=(20.5593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0493,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32672,0.0631472,-0.000109997,9.11124e-08,-2.91091e-11,76583,19.375], Tmin=(100,'K'), Tmax=(797.102,'K')), NASAPolynomial(coeffs=[11.7579,0.00918448,-5.4056e-06,1.09045e-09,-7.66069e-14,74971.4,-28.2615], Tmin=(797.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(635.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CCJ=O) + radical(CJ3)"""),
)

species(
    label = '[CH]=C1OC1[C]=O(3105)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u1 p0 c0 {2,D} {3,S}
6 C u1 p0 c0 {4,D} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (308.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95793,0.02814,2.06882e-05,-6.46403e-08,3.17493e-11,37135,18.6669], Tmin=(100,'K'), Tmax=(915.875,'K')), NASAPolynomial(coeffs=[19.9945,-0.00659995,5.46848e-06,-1.06847e-09,6.77485e-14,31984.3,-76.8461], Tmin=(915.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'O=[C]C=C1[CH]O1(3154)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u1 p0 c0 {1,S} {3,S} {7,S}
5 C u0 p0 c0 {3,D} {6,S} {8,S}
6 C u1 p0 c0 {2,D} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (178.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72713,0.0349558,2.32106e-06,-4.18692e-08,2.17799e-11,21589.8,16.8296], Tmin=(100,'K'), Tmax=(951.488,'K')), NASAPolynomial(coeffs=[19.4161,-0.00206219,1.8049e-06,-2.56879e-10,7.92611e-15,16533.2,-76.5133], Tmin=(951.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(methyleneoxirane) + radical(C=CCJO) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=[C]C1[CH]C1=O(3131)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {1,D} {3,S} {4,S}
6 C u1 p0 c0 {2,D} {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (248.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75955,0.0147234,3.65234e-05,-6.26342e-08,2.63037e-11,29925.5,20.3975], Tmin=(100,'K'), Tmax=(939.767,'K')), NASAPolynomial(coeffs=[12.2008,0.00643176,-1.14858e-06,2.03097e-10,-1.93727e-14,26742.7,-32.0575], Tmin=(939.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]=C=C[C]=O(525)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,D} {4,S} {6,S}
3 C u0 p0 c0 {2,D} {5,D}
4 C u1 p0 c0 {1,D} {2,S}
5 C u1 p0 c0 {3,D} {7,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (368.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1855,455,950,3120,650,792.5,1650,180.286,181.818,187.135,1693.12],'cm^-1')),
        HinderedRotor(inertia=(0.00829392,'amu*angstrom^2'), symmetry=1, barrier=(2.91538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0579,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5816,0.0342376,-4.28922e-05,3.11336e-08,-9.52278e-12,44355.1,15.702], Tmin=(100,'K'), Tmax=(783.921,'K')), NASAPolynomial(coeffs=[6.13377,0.0161129,-8.21218e-06,1.64147e-09,-1.17679e-13,43798.2,-0.570477], Tmin=(783.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C([O])C=C=O(1607)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p0 c0 {3,S} {6,D} {7,S}
5 C u1 p0 c0 {3,D} {8,S}
6 C u0 p0 c0 {2,D} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (191.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.76064,'amu*angstrom^2'), symmetry=1, barrier=(17.4886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69693,0.0510066,-7.12192e-05,4.97097e-08,-1.34031e-11,23123.1,19.2487], Tmin=(100,'K'), Tmax=(918.811,'K')), NASAPolynomial(coeffs=[11.2242,0.00953169,-3.51218e-06,5.85078e-10,-3.72633e-14,21372.3,-25.9085], Tmin=(918.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'O=C1[CH]C(=O)[CH]1(3155)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 O u0 p2 c0 {6,D}
3 C u1 p0 c0 {5,S} {6,S} {7,S}
4 C u1 p0 c0 {5,S} {6,S} {8,S}
5 C u0 p0 c0 {1,D} {3,S} {4,S}
6 C u0 p0 c0 {2,D} {3,S} {4,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (138.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26611,0.0270568,1.04873e-05,-3.49837e-08,1.5251e-11,16689.5,13.8289], Tmin=(100,'K'), Tmax=(1007.69,'K')), NASAPolynomial(coeffs=[12.8172,0.0113342,-5.04893e-06,1.05677e-09,-8.1686e-14,13234.8,-43.7452], Tmin=(1007.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(C=OCJC=O) + radical(C=OCJC=O)"""),
)

species(
    label = '[CH2]C([O])=[C][C]=O(1588)',
    structure = adjacencyList("""multiplicity 5
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 C u1 p0 c0 {3,D} {6,S}
6 C u1 p0 c0 {2,D} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (369.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.680956,'amu*angstrom^2'), symmetry=1, barrier=(15.6565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679327,'amu*angstrom^2'), symmetry=1, barrier=(15.6191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15515,0.0670563,-0.000114986,9.48724e-08,-3.01341e-11,44524.5,22.9319], Tmin=(100,'K'), Tmax=(829.761,'K')), NASAPolynomial(coeffs=[11.9655,0.0108005,-5.80098e-06,1.13139e-09,-7.79287e-14,42873.1,-26.3451], Tmin=(829.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C([O])=C=[C]O(3156)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {5,S} {8,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {2,S} {4,D} {6,S}
4 C u0 p0 c0 {3,D} {5,D}
5 C u1 p0 c0 {1,S} {4,D}
6 C u2 p0 c0 {3,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (471.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,540,610,2055,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10323,'amu*angstrom^2'), symmetry=1, barrier=(48.3575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10344,'amu*angstrom^2'), symmetry=1, barrier=(48.3623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922505,0.0595293,-7.07655e-05,4.04611e-08,-8.66271e-12,56836.4,24.7206], Tmin=(100,'K'), Tmax=(1293.94,'K')), NASAPolynomial(coeffs=[15.8381,0.00615172,-4.61535e-07,-1.023e-10,1.32328e-14,53584.9,-48.731], Tmin=(1293.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1[CH]C(=O)O1(3157)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u1 p0 c0 {1,S} {4,S} {6,S}
4 C u1 p0 c0 {3,S} {5,S} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (515.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75763,0.0182373,2.35514e-05,-4.83545e-08,2.18548e-11,62047,22.1705], Tmin=(100,'K'), Tmax=(897.813,'K')), NASAPolynomial(coeffs=[10.3203,0.00878575,-1.15954e-06,6.90686e-11,-3.29628e-15,59712,-18.9409], Tmin=(897.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C2CsJOC(O)) + radical(CCJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])[C][C]=O(3158)',
    structure = adjacencyList("""multiplicity 7
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u2 p0 c0 {3,S} {6,S}
5 C u2 p0 c0 {3,S} {8,S}
6 C u1 p0 c0 {2,D} {4,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (887.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1855,455,950,383.858,383.858,383.858,383.858,383.858,383.858],'cm^-1')),
        HinderedRotor(inertia=(0.0106578,'amu*angstrom^2'), symmetry=1, barrier=(1.11439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106619,'amu*angstrom^2'), symmetry=1, barrier=(1.11482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106621,'amu*angstrom^2'), symmetry=1, barrier=(1.11483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40772,0.0615181,-0.00010254,8.52562e-08,-2.72644e-11,106792,24.8353], Tmin=(100,'K'), Tmax=(867.642,'K')), NASAPolynomial(coeffs=[10.0835,0.0137653,-6.5752e-06,1.21811e-09,-8.13868e-14,105578,-14.1065], Tmin=(867.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(887.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C1=C[C]OO1(3159)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {5,S}
3 C u0 p0 c0 {1,S} {4,D} {6,S}
4 C u0 p0 c0 {3,D} {5,S} {7,S}
5 C u2 p0 c0 {2,S} {4,S}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (795.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13808,0.0330235,-6.76048e-06,-1.41018e-08,7.19354e-12,95754.3,20.0349], Tmin=(100,'K'), Tmax=(1052.17,'K')), NASAPolynomial(coeffs=[10.8916,0.0165626,-7.26841e-06,1.411e-09,-1.01808e-13,92981.4,-27.0645], Tmin=(1052.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(795.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxolene) + radical(CH2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])C=[C]O1(3122)',
    structure = adjacencyList("""multiplicity 5
1 O u0 p2 c0 {3,S} {5,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u1 p0 c0 {1,S} {4,D}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (670.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49548,0.0356059,1.70569e-05,-7.36363e-08,3.86196e-11,80756.3,20.0325], Tmin=(100,'K'), Tmax=(895.735,'K')), NASAPolynomial(coeffs=[24.3231,-0.0131892,9.77406e-06,-1.97973e-09,1.336e-13,74534.8,-99.4847], Tmin=(895.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CC(C)(O)OJ) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C][O](3112)',
    structure = adjacencyList("""multiplicity 6
1 O u1 p2 c0 {3,S}
2 C u2 p0 c0 {3,S} {4,S}
3 C u2 p0 c0 {1,S} {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (886.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,1121.29,1121.3,1121.31,1121.31],'cm^-1')),
        HinderedRotor(inertia=(0.217864,'amu*angstrom^2'), symmetry=1, barrier=(5.00912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88142,0.00967607,-6.02809e-06,1.41549e-09,-1.16526e-13,106618,9.68765], Tmin=(100,'K'), Tmax=(3543.44,'K')), NASAPolynomial(coeffs=[2.07116,0.00784258,-3.61075e-06,6.51916e-10,-4.08687e-14,110335,24.1457], Tmin=(3543.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(886.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CCOJ) + radical(CCJ2_triplet) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C([O])=[C][CH][O](3160)',
    structure = adjacencyList("""multiplicity 7
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {5,D} {6,S}
4 C u1 p0 c0 {2,S} {5,S} {7,S}
5 C u1 p0 c0 {3,D} {4,S}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (691.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,1685,370,495.696,496.975,498.536,501.047,501.459,503.888],'cm^-1')),
        HinderedRotor(inertia=(0.000669433,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290075,'amu*angstrom^2'), symmetry=1, barrier=(51.7513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12248,0.0425601,-3.93423e-05,2.06696e-08,-4.62144e-12,83263.1,24.0695], Tmin=(100,'K'), Tmax=(1048.53,'K')), NASAPolynomial(coeffs=[7.55668,0.0218294,-9.68558e-06,1.8136e-09,-1.25621e-13,82123.5,-2.40489], Tmin=(1048.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([O])[C]=C=O(3161)',
    structure = adjacencyList("""multiplicity 5
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u1 p0 c0 {3,S} {6,D}
5 C u2 p0 c0 {3,S} {8,S}
6 C u0 p0 c0 {2,D} {4,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (598.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,2120,512.5,787.5,323.006,323.009,323.015,323.015,323.017],'cm^-1')),
        HinderedRotor(inertia=(0.866787,'amu*angstrom^2'), symmetry=1, barrier=(64.1814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00500907,'amu*angstrom^2'), symmetry=1, barrier=(45.1917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32436,0.0635513,-0.000106078,8.80106e-08,-2.79107e-11,72053.9,21.2616], Tmin=(100,'K'), Tmax=(888.461,'K')), NASAPolynomial(coeffs=[10.2077,0.0140341,-6.40021e-06,1.15193e-09,-7.52752e-14,70851.2,-18.4295], Tmin=(888.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CC(C)OJ) + radical(CCCJ=C=O) + radical(CCJ2_triplet)"""),
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
    E0 = (-57.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (499.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (613.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-49.4367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (5.67917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (361.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (251.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (69.3435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (184.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (161.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (398.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (141.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (241.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (453.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (452.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-54.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-54.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-54.7922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (196.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (216.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-57.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-49.4367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (158.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (269.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (184.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (515.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (400.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (275.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (658.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (319.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (381.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['HCCO(20)', 'HCCO(20)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(7)', '[CH][C]=C[C]=O(519)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(179)', '[CH]C(=[CH])[O](3145)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]=C1[CH]C(=O)O1(3126)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]C(O)=C=C=O(3146)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]=C1OC1[C][O](3115)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(418.743,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 416.2 to 418.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]C1([O])C=C1[O](3147)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(308.738,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 305.8 to 308.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CO(14)', '[CH]C(=[CH])[O](3145)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(302000,'m^3/(mol*s)'), n=-5.25699e-09, Ea=(13.726,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_N-Sp-4R!H-3C',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_N-Sp-4R!H-3C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH](3148)', 'O=[C]C=C=O(470)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(16.08,'m^3/(mol*s)'), n=1.68, Ea=(2.63946,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Sp-2R!H=1R!H_N-2R!H->S_4R!H-u0_N-4R!H-inRing_N-2CO->C',), comment="""Estimated from node Root_N-3R-inRing_3R->C_Ext-1R!H-R_N-Sp-2R!H-=1R!H_Sp-2R!H=1R!H_N-2R!H->S_4R!H-u0_N-4R!H-inRing_N-2CO->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', '[CH]C([O])=C=C=O(3149)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.7698e-05,'m^3/(mol*s)'), n=3.50443, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.1710730160484305, var=0.169866296243597, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(6)', '[CH]C([O])=[C][C]=O(3150)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C([O])=C=C[O](3151)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(O)=[C][C]=O(3152)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH][C]=O(474)', '[CH][C]=O(474)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(6)', '[C]C([O])=C[C]=O(3153)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]=C1OC1[C]=O(3105)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['O=[C]C=C1[CH]O1(3154)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['O=[C]C1[CH]C1=O(3131)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HCCO(20)', '[CH][C]=O(474)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O(7)', '[CH]=C=C[C]=O(525)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.97583,'m^3/(mol*s)'), n=1.98069, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.0665490741671073, var=3.1459699634659253, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R_N-3BrOS-u1',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R_N-3BrOS-u1"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]=C([O])C=C=O(1607)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['O=C1[CH]C(=O)[CH]1(3155)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])=[C][C]=O(1588)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.97864e+06,'s^-1'), n=1.84533, Ea=(184.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SD;Y_rad_out;Cd_H_out_singleH] + [R3H_SD;Cd_rad_out;Cd_H_out_single] for rate rule [R3H_SD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]C([O])=C=[C]O(3156)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH][C]1[CH]C(=O)O1(3157)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra] for rate rule [R4_linear;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C([O])[C][C]=O(3158)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]C1=C[C]OO1(3159)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.85053e+09,'s^-1'), n=0.454312, Ea=(458.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 454.5 to 458.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C([O])=C[C]=O(1578)'],
    products = ['[CH]C1([O])C=[C]O1(3122)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(333.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;multiplebond_intra;radadd_intra_O] for rate rule [R5;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 330.1 to 333.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['HCCO(20)', '[CH][C][O](3112)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C([O])=[C][CH][O](3160)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C([O])[C]=C=O(3161)'],
    products = ['[CH]C([O])=C[C]=O(1578)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #887',
    isomers = [
        '[CH]C([O])=C[C]=O(1578)',
    ],
    reactants = [
        ('HCCO(20)', 'HCCO(20)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #887',
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

