species(
    label = '[CH]=CC(=O)[C]=O(1360)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,D} {3,S} {6,S}
5 C u1 p0 c0 {3,D} {8,S}
6 C u1 p0 c0 {2,D} {4,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (224.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,3120,650,792.5,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.522705,'amu*angstrom^2'), symmetry=1, barrier=(12.018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11681,'amu*angstrom^2'), symmetry=1, barrier=(25.6778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17022,0.0443142,-5.9188e-05,4.20312e-08,-1.11204e-11,27062.8,20.1803], Tmin=(100,'K'), Tmax=(638.017,'K')), NASAPolynomial(coeffs=[7.16348,0.0179797,-8.96004e-06,1.75796e-09,-1.24153e-13,26324.5,-2.45812], Tmin=(638.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C=C=O(813)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (63.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2110,2130,495,530,650,925,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2387.04,'J/mol'), sigma=(4.99307,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=372.85 K, Pc=43.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5355,0.024033,-4.26139e-05,3.87325e-08,-1.34306e-11,7697.1,25.5769], Tmin=(100,'K'), Tmax=(857.394,'K')), NASAPolynomial(coeffs=[4.78324,0.00773532,-3.93446e-06,7.52107e-10,-5.12483e-14,7525.26,16.3243], Tmin=(857.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.731,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""OCCO(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[C]=O(117)',
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
    label = '[CH]C=C=O(302)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u2 p0 c0 {2,S} {6,S}
4 C u0 p0 c0 {1,D} {2,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (292.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2120,512.5,787.5,243.591,248.35,249.237],'cm^-1')),
        HinderedRotor(inertia=(1.16574,'amu*angstrom^2'), symmetry=1, barrier=(50.5613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00496,0.0183622,-5.86376e-07,-8.19938e-09,3.38636e-12,35273.9,13.1075], Tmin=(100,'K'), Tmax=(1106.37,'K')), NASAPolynomial(coeffs=[5.9629,0.0150316,-6.05408e-06,1.11097e-09,-7.67694e-14,34168.7,-3.49851], Tmin=(1106.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=C1C=CC1=O(1372)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,D} {3,S} {6,S}
5 C u0 p0 c0 {3,D} {6,S} {8,S}
6 C u0 p0 c0 {2,D} {4,S} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-61.6105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74105,0.0217604,3.46192e-06,-1.58338e-08,6.12344e-12,-7359.79,14.1742], Tmin=(100,'K'), Tmax=(1127.7,'K')), NASAPolynomial(coeffs=[8.09217,0.0164637,-7.69452e-06,1.5221e-09,-1.09705e-13,-9436.77,-16.1428], Tmin=(1127.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.6105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C[C]1OC1=O(1373)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 C u0 p0 c0 {3,S} {6,D} {7,S}
6 C u1 p0 c0 {5,D} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (293.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42019,0.0388739,-5.47722e-05,5.02207e-08,-1.89428e-11,35359.1,19.5116], Tmin=(100,'K'), Tmax=(791.414,'K')), NASAPolynomial(coeffs=[4.10725,0.0231656,-1.13881e-05,2.20911e-09,-1.54401e-13,35317,13.188], Tmin=(791.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(2(co)oxirane) + radical(C2CsJO) + radical(Cds_P)"""),
)

species(
    label = 'O=C=C1[CH][CH]O1(1370)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {6,D}
4 C u1 p0 c0 {3,S} {5,S} {7,S}
5 C u1 p0 c0 {1,S} {4,S} {8,S}
6 C u0 p0 c0 {2,D} {3,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (235.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51613,0.0361026,1.44539e-05,-6.44547e-08,3.26974e-11,28422.3,14.4294], Tmin=(100,'K'), Tmax=(921.332,'K')), NASAPolynomial(coeffs=[22.4531,-0.00654906,5.34441e-06,-1.02542e-09,6.33619e-14,22516.6,-95.9765], Tmin=(921.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + ring(Cyclobutane) + radical(C=CCJCO) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C1([C]=O)C=C1(1374)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {3,S} {4,D} {8,S}
6 C u1 p0 c0 {2,D} {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (374.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98532,0.0498397,-8.34761e-05,7.53132e-08,-2.62703e-11,45071.1,19.7033], Tmin=(100,'K'), Tmax=(849.21,'K')), NASAPolynomial(coeffs=[6.23318,0.0188184,-9.22933e-06,1.75536e-09,-1.19909e-13,44746.7,2.24216], Tmin=(849.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CC(C)(C=O)OJ) + radical(CC(C)(O)CJ=O)"""),
)

species(
    label = '[CH]=C=C(O)[C]=O(1375)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {7,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,D}
5 C u1 p0 c0 {2,D} {3,S}
6 C u1 p0 c0 {4,D} {8,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (146.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,540,610,2055,1855,455,950,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.19055,'amu*angstrom^2'), symmetry=1, barrier=(27.3732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21782,0.0576063,-7.60577e-05,4.65772e-08,-1.07836e-11,17740.4,19.3548], Tmin=(100,'K'), Tmax=(1077.06,'K')), NASAPolynomial(coeffs=[15.7154,0.00376518,-1.07431e-06,1.64872e-10,-1.0748e-14,14617.4,-51.6637], Tmin=(1077.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O) + radical(C=C=CJ)"""),
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
    label = '[O]C#C[O](814)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
"""),
    E0 = (9.55717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,344.097,344.099,344.1],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09811,0.0228497,-4.17927e-05,3.86217e-08,-1.34608e-11,1179.06,7.869], Tmin=(100,'K'), Tmax=(871.084,'K')), NASAPolynomial(coeffs=[4.9831,0.00749284,-3.80902e-06,7.20359e-10,-4.85535e-14,1104.9,0.494495], Tmin=(871.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.55717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""OCCO""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'C#CC(=O)[C]=O(1376)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,D} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {6,T}
5 C u1 p0 c0 {2,D} {3,S}
6 C u0 p0 c0 {4,T} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (143.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2175,525,1855,455,950,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.76982,'amu*angstrom^2'), symmetry=1, barrier=(40.6917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77326,'amu*angstrom^2'), symmetry=1, barrier=(40.7706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0493,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0436,0.0450846,-6.76733e-05,5.24468e-08,-1.57881e-11,17271.5,17.9675], Tmin=(100,'K'), Tmax=(894.079,'K')), NASAPolynomial(coeffs=[8.78505,0.0112176,-4.63624e-06,8.06798e-10,-5.21733e-14,16214.1,-12.9728], Tmin=(894.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[CH](149)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,D} {3,S}
2 C u1 p0 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1227.31,2941.71],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88058,-0.000843434,2.04094e-05,-2.51944e-08,9.47575e-12,71059.3,4.72115], Tmin=(100,'K'), Tmax=(935.367,'K')), NASAPolynomial(coeffs=[4.92141,0.00358272,-9.24433e-07,1.57271e-10,-1.19539e-14,70476.2,-2.30656], Tmin=(935.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]C([O])=C=O(1377)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {5,S} {6,D}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u1 p0 c0 {3,S} {4,D}
6 C u0 p0 c0 {2,D} {3,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (204.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58826,0.0564729,-9.05044e-05,7.22954e-08,-2.17457e-11,24692.5,18.8525], Tmin=(100,'K'), Tmax=(961.111,'K')), NASAPolynomial(coeffs=[9.64202,0.0120589,-4.18324e-06,6.24364e-10,-3.48385e-14,23647.6,-17.0648], Tmin=(961.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsOs) + group(Cd-Cd(CCO)H) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC([O])=C[O](1378)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {7,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (127.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605977,0.054746,-6.14816e-05,3.08791e-08,-5.58328e-12,15448.7,21.0975], Tmin=(100,'K'), Tmax=(1636.21,'K')), NASAPolynomial(coeffs=[17.7522,-0.00135756,2.95723e-06,-6.75382e-10,4.76653e-14,11736.8,-64.2623], Tmin=(1636.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
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
    label = '[CH]C=C=C=O(1379)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 C u0 p0 c0 {3,D} {4,S} {6,S}
3 C u0 p0 c0 {2,D} {5,D}
4 C u2 p0 c0 {2,S} {7,S}
5 C u0 p0 c0 {1,D} {3,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (421.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,540,610,2055,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15637,'amu*angstrom^2'), symmetry=1, barrier=(49.5791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0579,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50245,0.0299006,-1.87648e-05,5.85794e-09,-7.51457e-13,50724.8,16.1895], Tmin=(100,'K'), Tmax=(1750.16,'K')), NASAPolynomial(coeffs=[8.87316,0.0153403,-6.28576e-06,1.10446e-09,-7.24526e-14,48494.8,-18.1113], Tmin=(1750.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cdd-(Cdd-O2d)Cds) + missing(Cdd-CddO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=C=C1C=CO1(1362)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {6,D}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {1,S} {4,D} {8,S}
6 C u0 p0 c0 {2,D} {3,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (96.5354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42496,0.0271128,-3.29861e-06,-1.95379e-08,1.0682e-11,11674.1,14.6594], Tmin=(100,'K'), Tmax=(947.929,'K')), NASAPolynomial(coeffs=[11.4012,0.00782799,-2.20279e-06,3.82529e-10,-2.86144e-14,9137,-32.5718], Tmin=(947.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsOs) + group(Cd-Cd(CCO)H) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + ring(Oxetene)"""),
)

species(
    label = 'O=C1[CH][CH]C1=O(1380)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 O u0 p2 c0 {6,D}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {1,D} {4,S} {6,S}
6 C u0 p0 c0 {2,D} {3,S} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (205.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84024,0.0182007,1.39686e-05,-2.75597e-08,1.04757e-11,24782.4,19.2218], Tmin=(100,'K'), Tmax=(1056.48,'K')), NASAPolynomial(coeffs=[8.24496,0.0153256,-6.92075e-06,1.37953e-09,-1.01096e-13,22658.9,-11.7951], Tmin=(1056.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclobutane) + radical(CCJCC=O) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]=CC1=[C]OO1(1381)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p0 c0 {3,S} {6,D} {7,S}
5 C u1 p0 c0 {2,S} {3,D}
6 C u1 p0 c0 {4,D} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (677.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12186,0.0350377,-2.31053e-05,1.53103e-09,2.60057e-12,81548.5,20.5815], Tmin=(100,'K'), Tmax=(1027.66,'K')), NASAPolynomial(coeffs=[11.8546,0.00950713,-3.87065e-06,7.49887e-10,-5.49155e-14,78895.8,-29.8126], Tmin=(1027.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[O]C1=[C]OC=C1(1369)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {6,S}
2 O u1 p2 c0 {5,S}
3 C u0 p0 c0 {4,D} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {8,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {1,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (116.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75455,0.0189867,1.48904e-05,-3.35434e-08,1.42107e-11,14059,18.3353], Tmin=(100,'K'), Tmax=(977.024,'K')), NASAPolynomial(coeffs=[9.77003,0.0107514,-3.91862e-06,7.52285e-10,-5.63484e-14,11710.3,-20.3512], Tmin=(977.024,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Furan) + radical(C=C(C)OJ) + radical(C=CJO)"""),
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
    E0 = (52.7479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (503.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (4.01407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (164.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (133.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (145.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (179.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (9.1524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (58.7912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (140.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (425.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (371.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (163.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (119.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (435.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (4.01407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (133.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (448.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (29.1023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['O=C=C=O(813)', 'C2H2(22)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(57.0182,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 57.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=O(117)', '[CH]C=C=O(302)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['O=C1C=CC1=O(1372)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['[CH]=C[C]1OC1=O(1373)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(169.231,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_De;radadd_intra_CO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['O=C=C1[CH][CH]O1(1370)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['[O]C1([C]=O)C=C1(1374)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(150.016,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=C(O)[C]=O(1375)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(261.378,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CO(14)', '[CH]C=C=O(302)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(86.1,'m^3/(mol*s)'), n=1.36, Ea=(63.6986,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-3C-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-3C-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C#C[O](814)', 'C2H2(22)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.63031e-05,'m^3/(mol*s)'), n=3.04336, Ea=(60.2171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', 'C#CC(=O)[C]=O(1376)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1348.96,'m^3/(mol*s)'), n=1.432, Ea=(14.6518,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C=C=O(813)', '[CH]=[CH](149)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(647.032,'m^3/(mol*s)'), n=1.13767, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.011839139166081196, var=0.8671401158115553, Tref=1000.0, N=160, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C#C[O](814)', '[CH]=[CH](149)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.12792e+08,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['C=[C]C([O])=C=O(1377)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['C#CC([O])=C[O](1378)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(17055.2,'s^-1'), n=2.44296, Ea=(123.676,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;XH_out] for rate rule [R3H_SS;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O(7)', '[CH]C=C=C=O(1379)'],
    products = ['[CH]=CC(=O)[C]=O(1360)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['O=C=C1C=CO1(1362)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['O=C1[CH][CH]C1=O(1380)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['[CH]=CC1=[C]OO1(1381)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.8957e+11,'s^-1'), n=0.254, Ea=(452.928,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 451.5 to 452.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC(=O)[C]=O(1360)'],
    products = ['[O]C1=[C]OC=C1(1369)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.17816e+10,'s^-1'), n=0.292763, Ea=(33.3726,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_MS;multiplebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #711',
    isomers = [
        '[CH]=CC(=O)[C]=O(1360)',
    ],
    reactants = [
        ('O=C=C=O(813)', 'C2H2(22)'),
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

