species(
    label = 'O=[C]OC(=O)[C]=O(1308)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {7,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,D} {6,S}
6 C u1 p0 c0 {3,D} {5,S}
7 C u1 p0 c0 {1,S} {4,D}
"""),
    E0 = (-223.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1850,1860,440,470,900,1000,261.681,261.681,261.682,261.682,847.49,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00246181,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314065,'amu*angstrom^2'), symmetry=1, barrier=(15.2613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622248,'amu*angstrom^2'), symmetry=1, barrier=(30.2369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56138,0.0592822,-0.000107158,9.41784e-08,-3.17281e-11,-26834.3,23.7332], Tmin=(100,'K'), Tmax=(836.034,'K')), NASAPolynomial(coeffs=[9.41003,0.0126355,-7.147e-06,1.41602e-09,-9.8268e-14,-27828.8,-10.8255], Tmin=(836.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + group(Cds-OdOsH) + radical(OC=OCJ=O) + radical((O)CJOC)"""),
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
    label = 'O=[C]O[C]=O(290)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,D}
5 C u1 p0 c0 {1,S} {3,D}
"""),
    E0 = (-91.7908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1850,1860,440,470,900,1000,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.366629,'amu*angstrom^2'), symmetry=1, barrier=(12.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364921,'amu*angstrom^2'), symmetry=1, barrier=(12.334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0195,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.87,'J/mol'), sigma=(5.46152,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.17 K, Pc=48.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47246,0.037474,-6.98113e-05,6.12863e-08,-2.02076e-11,-10988.6,15.1735], Tmin=(100,'K'), Tmax=(882.864,'K')), NASAPolynomial(coeffs=[7.43584,0.00673308,-3.55935e-06,6.69641e-10,-4.43946e-14,-11543.3,-6.33194], Tmin=(882.864,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.7908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(99.7737,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-OdOsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical((O)CJOC)"""),
)

species(
    label = '[O]C([O])=C=O(295)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
"""),
    E0 = (-100.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2120,512.5,787.5,180,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0195,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3750.28,'J/mol'), sigma=(6.30935,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=585.79 K, Pc=33.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16765,0.0177643,4.73257e-05,-1.02219e-07,4.85607e-11,-12010.1,14.1784], Tmin=(100,'K'), Tmax=(899.698,'K')), NASAPolynomial(coeffs=[24.743,-0.0222738,1.34934e-05,-2.61717e-09,1.73914e-13,-18514.1,-105.918], Tmin=(899.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + missing(Cdd-CdO2d) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=O(213)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3422.77,'J/mol'), sigma=(5.42046,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=534.63 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186753,1.03202e-05,-1.52373e-08,5.80537e-12,3804.5,8.40407], Tmin=(100,'K'), Tmax=(1021.26,'K')), NASAPolynomial(coeffs=[6.36178,0.000422737,-4.06602e-07,1.525e-10,-1.51974e-14,2816.75,-6.43921], Tmin=(1021.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = 'O=C1OC(=O)C1=O(1311)',
    structure = adjacencyList("""1 O u0 p2 c0 {6,S} {7,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,D} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 C u0 p0 c0 {1,S} {4,D} {5,S}
"""),
    E0 = (-435.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50814,0.0348088,-3.60791e-05,1.85737e-08,-3.88032e-12,-52341,18.0335], Tmin=(100,'K'), Tmax=(1135.73,'K')), NASAPolynomial(coeffs=[8.7909,0.0126809,-6.85369e-06,1.41834e-09,-1.03991e-13,-53768,-13.0768], Tmin=(1135.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)O2s) + ring(Cyclobutane)"""),
)

species(
    label = 'O=[C]O[C]1OC1=O(1318)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 C u1 p0 c0 {2,S} {4,D}
"""),
    E0 = (-150.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24264,0.0435301,-6.74193e-05,6.35605e-08,-2.44317e-11,-18039.3,19.9786], Tmin=(100,'K'), Tmax=(765.509,'K')), NASAPolynomial(coeffs=[4.64903,0.0232192,-1.24604e-05,2.4953e-09,-1.77339e-13,-18181,10.4929], Tmin=(765.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + group(Cds-OdOsH) + ring(2(co)oxirane) + radical(Cs_P) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C][C]1OC(=O)O1(1319)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {2,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {3,D}
7 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-220.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84923,0.00701448,6.50145e-05,-9.3272e-08,3.61739e-11,-26480.7,20.8634], Tmin=(100,'K'), Tmax=(979.998,'K')), NASAPolynomial(coeffs=[14.7181,0.00478518,-2.31112e-06,6.48807e-10,-6.12862e-14,-31026.2,-47.4801], Tmin=(979.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(Cds-OdOsOs) + ring(Cyclobutane) + radical(Cs_P) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C1([C]=O)OC1=O(1320)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-126.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13094,0.0562249,-7.41845e-05,4.41995e-08,-9.73459e-12,-15132.8,21.112], Tmin=(100,'K'), Tmax=(1244.97,'K')), NASAPolynomial(coeffs=[16.8104,-0.000400531,1.56885e-06,-3.96927e-10,3.0256e-14,-18552.7,-56.0233], Tmin=(1244.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cds-OdCsOs) + group(Cds-OdCsH) + ring(2(co)oxirane) + radical(C=OCOJ) + radical(C=OCCJ=O)"""),
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
    label = 'O=[C]O[C]=C=O(1321)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u1 p0 c0 {1,S} {6,D}
5 C u1 p0 c0 {1,S} {2,D}
6 C u0 p0 c0 {3,D} {4,D}
"""),
    E0 = (114.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,1855,455,950,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48681,'amu*angstrom^2'), symmetry=1, barrier=(34.1847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998989,'amu*angstrom^2'), symmetry=1, barrier=(22.9687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0301,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95562,0.0491495,-8.63505e-05,7.39874e-08,-2.42957e-11,13862,18.8841], Tmin=(100,'K'), Tmax=(851.231,'K')), NASAPolynomial(coeffs=[8.79472,0.0101044,-5.3749e-06,1.03604e-09,-7.06852e-14,12947.9,-11.5391], Tmin=(851.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(Cds-OdOsH) + missing(Cdd-CdO2d) + radical(C=CJO) + radical((O)CJOC)"""),
)

species(
    label = 'O=C=C1OC(=O)O1(1322)',
    structure = adjacencyList("""1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {7,D}
6 C u0 p0 c0 {1,S} {2,S} {3,D}
7 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (-379.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53109,0.022023,7.61867e-05,-1.46157e-07,6.57362e-11,-45498.7,15.5167], Tmin=(100,'K'), Tmax=(922.651,'K')), NASAPolynomial(coeffs=[31.238,-0.0230564,1.33827e-05,-2.44342e-09,1.51441e-13,-54543.6,-144.719], Tmin=(922.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + group(Cds-OdOsOs) + missing(Cdd-CdO2d) + ring(Cyclobutane)"""),
)

species(
    label = '[O][C]1C(=O)OC1=O(1323)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {6,S} {7,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {2,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 C u0 p0 c0 {1,S} {4,D} {5,S}
"""),
    E0 = (-166.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61048,0.0272016,-2.01816e-05,7.20068e-09,-1.0221e-12,-19945.1,22.5627], Tmin=(100,'K'), Tmax=(1651.64,'K')), NASAPolynomial(coeffs=[9.79919,0.00979168,-4.37013e-06,8.18531e-10,-5.60662e-14,-22319.8,-15.7257], Tmin=(1651.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsOs) + ring(Cyclobutane) + radical(C=OCOJ) + radical(C2CsJOH)"""),
)

species(
    label = 'O=[C]OC1=[C]OO1(1324)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u1 p0 c0 {3,S} {5,D}
7 C u1 p0 c0 {2,S} {4,D}
"""),
    E0 = (295.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07089,0.0482103,-6.95404e-05,4.89537e-08,-1.11647e-11,35593.9,20.6921], Tmin=(100,'K'), Tmax=(577.24,'K')), NASAPolynomial(coeffs=[7.40279,0.0183824,-1.05314e-05,2.1702e-09,-1.57029e-13,34859.7,-3.12876], Tmin=(577.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdOsH) + ring(Cyclobutene) + radical(C=CJO) + radical((O)CJOC)"""),
)

species(
    label = '[O]C1=[C]OC(=O)O1(1325)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 C u0 p0 c0 {1,S} {3,S} {7,D}
7 C u1 p0 c0 {2,S} {6,D}
"""),
    E0 = (-228.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16293,0.0236176,2.76193e-05,-6.52473e-08,2.9534e-11,-27450,18.7601], Tmin=(100,'K'), Tmax=(948.006,'K')), NASAPolynomial(coeffs=[18.3763,-0.00253676,2.14165e-06,-3.12165e-10,1.05544e-14,-32422.9,-68.6094], Tmin=(948.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + ring(Cyclopentane) + radical(C=COJ) + radical(C=CJO)"""),
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
    E0 = (-160.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (401.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (410.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-152.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (8.45383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (81.3924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-63.7097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-120.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-121.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-153.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (158.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (104.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (420.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-152.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-34.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (358.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-105.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['CO2(15)', 'O=C=C=O(813)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=O(117)', '[O]C([O])=C=O(295)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(117)', 'O=[C]O[C]=O(290)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['O=C1OC(=O)C1=O(1311)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['O=[C]O[C]1OC1=O(1318)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(169.231,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_CO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['O=[C][C]1OC(=O)O1(1319)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['[O]C1([C]=O)OC1=O(1320)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.1935e+10,'s^-1'), n=0.672833, Ea=(97.0679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_CO]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Exocyclic
Ea raised from 96.8 to 97.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CO(14)', 'O=[C]O[C]=O(290)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.659862,'m^3/(mol*s)'), n=2.07304, Ea=(26.8124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17748312369238067, var=0.6615994821919932, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CO(14)', '[O]C([O])=C=O(295)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(68.2,'m^3/(mol*s)'), n=8.73864e-09, Ea=(34.519,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CO2(15)', '[O]C#C[O](814)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.63031e-05,'m^3/(mol*s)'), n=3.04336, Ea=(177.025,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(213)', 'O=C=C=O(813)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.64086e-07,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(213)', '[O]C#C[O](814)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.95326e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O(7)', 'O=[C]O[C]=C=O(1321)'],
    products = ['O=[C]OC(=O)[C]=O(1308)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['O=C=C1OC(=O)O1(1322)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['[O][C]1C(=O)OC1=O(1323)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['O=[C]OC1=[C]OO1(1324)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.8957e+11,'s^-1'), n=0.254, Ea=(519.211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 518.3 to 519.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]OC(=O)[C]=O(1308)'],
    products = ['[O]C1=[C]OC(=O)O1(1325)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.85053e+09,'s^-1'), n=0.454312, Ea=(54.8171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;multiplebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #702',
    isomers = [
        'O=[C]OC(=O)[C]=O(1308)',
    ],
    reactants = [
        ('CO2(15)', 'O=C=C=O(813)'),
        ('CO(14)', 'O=[C]O[C]=O(290)'),
        ('CO(14)', '[O]C([O])=C=O(295)'),
        ('[O][C]=O(213)', 'O=C=C=O(813)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #702',
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

