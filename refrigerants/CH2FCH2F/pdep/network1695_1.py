species(
    label = '[CH]=C(O)OC#C(1479)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u1 p0 c0 {3,D} {8,S}
5 C u0 p0 c0 {1,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (269.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3120,650,792.5,1650,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.07975,'amu*angstrom^2'), symmetry=1, barrier=(24.8256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08073,'amu*angstrom^2'), symmetry=1, barrier=(24.848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08002,'amu*angstrom^2'), symmetry=1, barrier=(24.8318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10946,0.0638236,-8.98894e-05,6.10557e-08,-1.5915e-11,32469.1,19.0638], Tmin=(100,'K'), Tmax=(951.347,'K')), NASAPolynomial(coeffs=[14.1132,0.00915207,-3.69342e-06,6.56716e-10,-4.40059e-14,29994.8,-43.0238], Tmin=(951.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = 'HCCO(21)',
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
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[C]=C(O)OC#C(3930)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {6,D}
4 C u0 p0 c0 {1,S} {5,T}
5 C u0 p0 c0 {4,T} {8,S}
6 C u2 p0 c0 {3,D}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (580.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01589,'amu*angstrom^2'), symmetry=1, barrier=(23.3572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0155,'amu*angstrom^2'), symmetry=1, barrier=(23.3483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0161,'amu*angstrom^2'), symmetry=1, barrier=(23.362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11144,0.0671194,-0.000112114,8.95459e-08,-2.73603e-11,69871.2,19.0298], Tmin=(100,'K'), Tmax=(868.539,'K')), NASAPolynomial(coeffs=[12.6165,0.00949027,-4.56713e-06,8.40447e-10,-5.56071e-14,68047.8,-33.8456], Tmin=(868.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'OC1=CC=[C]O1(3931)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {3,S} {9,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {5,S} {7,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u1 p0 c0 {1,S} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (35.4089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53027,0.0260244,-2.72206e-06,-1.27343e-08,5.84419e-12,4316.97,18.1868], Tmin=(100,'K'), Tmax=(1078.63,'K')), NASAPolynomial(coeffs=[9.03356,0.0153879,-6.67685e-06,1.29661e-09,-9.33157e-14,2129.86,-17.3151], Tmin=(1078.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.4089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Furan) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1C=C(O)O1(3932)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {5,S} {7,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u1 p0 c0 {5,D} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (198.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31092,0.0478203,-4.9837e-05,2.4686e-08,-4.56891e-12,24006.3,16.9932], Tmin=(100,'K'), Tmax=(1503.82,'K')), NASAPolynomial(coeffs=[15.1765,0.00460973,-4.226e-07,-1.91265e-11,3.45572e-15,20551.7,-53.1776], Tmin=(1503.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Oxetene) + radical(Cds_P)"""),
)

species(
    label = 'C#CO(429)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {5,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,653.448,840.947],'cm^-1')),
        HinderedRotor(inertia=(0.2322,'amu*angstrom^2'), symmetry=1, barrier=(5.33873,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'OH(7)',
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
    label = 'C#COC#C(3933)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {4,T}
3 C u0 p0 c0 {1,S} {5,T}
4 C u0 p0 c0 {2,T} {6,S}
5 C u0 p0 c0 {3,T} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (388.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,180],'cm^-1')),
        HinderedRotor(inertia=(1.6047,'amu*angstrom^2'), symmetry=1, barrier=(36.8952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59463,'amu*angstrom^2'), symmetry=1, barrier=(36.6636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0579,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03948,0.0406272,-4.7466e-05,2.66876e-08,-5.78225e-12,46794.5,12.2972], Tmin=(100,'K'), Tmax=(1138.69,'K')), NASAPolynomial(coeffs=[11.6997,0.00669264,-2.76394e-06,5.16019e-10,-3.62717e-14,44594.5,-35.5626], Tmin=(1138.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Ct-CtOs) + group(Ct-CtOs) + group(Ct-CtH) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=[C]O(432)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.92,1059.79],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10435e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.311,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C2H(22)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]C(=O)O(3934)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {6,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (2.71927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,717.252,717.252,717.252,717.252,717.252,717.252,717.252],'cm^-1')),
        HinderedRotor(inertia=(0.00730754,'amu*angstrom^2'), symmetry=1, barrier=(2.66773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00730754,'amu*angstrom^2'), symmetry=1, barrier=(2.66773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.15932,0.0134425,6.99649e-06,-1.70171e-08,6.7659e-12,361.656,13.2556], Tmin=(100,'K'), Tmax=(1050.37,'K')), NASAPolynomial(coeffs=[7.43206,0.00914904,-3.97751e-06,8.04797e-10,-5.99604e-14,-1196.68,-10.7132], Tmin=(1050.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.71927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C]OC#C(3935)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 C u1 p0 c0 {1,S} {4,D}
3 C u0 p0 c0 {1,S} {5,T}
4 C u1 p0 c0 {2,D} {7,S}
5 C u0 p0 c0 {3,T} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (665.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2175,525,3120,650,792.5,1650,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.25357,'amu*angstrom^2'), symmetry=1, barrier=(28.8221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25441,'amu*angstrom^2'), symmetry=1, barrier=(28.8413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0579,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26086,0.039219,-5.10512e-05,3.39548e-08,-8.88437e-12,80160.8,17.5811], Tmin=(100,'K'), Tmax=(938.18,'K')), NASAPolynomial(coeffs=[9.09881,0.0100644,-4.43692e-06,8.30391e-10,-5.74527e-14,78877.8,-14.9716], Tmin=(938.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=O)OC#C(3936)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {5,S}
4 C u0 p0 c0 {1,S} {6,T}
5 C u2 p0 c0 {3,S} {7,S}
6 C u0 p0 c0 {4,T} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (327.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89537,0.0410592,-3.66548e-05,1.56485e-08,-2.6201e-12,39422,18.9167], Tmin=(100,'K'), Tmax=(1439.3,'K')), NASAPolynomial(coeffs=[12.8528,0.0106067,-4.9178e-06,9.48219e-10,-6.66981e-14,36267.9,-37.9368], Tmin=(1439.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]#COC(=[CH])O(3937)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u1 p0 c0 {3,D} {8,S}
5 C u0 p0 c0 {1,S} {6,T}
6 C u1 p0 c0 {5,T}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (606.243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3120,650,792.5,1650,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.947133,'amu*angstrom^2'), symmetry=1, barrier=(21.7765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.947205,'amu*angstrom^2'), symmetry=1, barrier=(21.7781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.947153,'amu*angstrom^2'), symmetry=1, barrier=(21.7769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07143,0.0681239,-0.000115294,9.18712e-08,-2.74985e-11,73016.2,19.9398], Tmin=(100,'K'), Tmax=(938.234,'K')), NASAPolynomial(coeffs=[12.6834,0.00830576,-3.17338e-06,4.88228e-10,-2.7298e-14,71291.2,-32.9222], Tmin=(938.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_P) + radical(Acetyl)"""),
)

species(
    label = 'C#COC([CH2])=O(1013)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (90.4817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180,180,180,756.278,756.294,756.374],'cm^-1')),
        HinderedRotor(inertia=(0.0965083,'amu*angstrom^2'), symmetry=1, barrier=(39.1529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0963241,'amu*angstrom^2'), symmetry=1, barrier=(39.1557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964423,'amu*angstrom^2'), symmetry=1, barrier=(39.1522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78591,0.0433021,-3.72489e-05,1.55609e-08,-2.56599e-12,10966.6,19.6326], Tmin=(100,'K'), Tmax=(1454.25,'K')), NASAPolynomial(coeffs=[12.9295,0.0126505,-5.6324e-06,1.06675e-09,-7.42583e-14,7725.53,-38.3015], Tmin=(1454.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.4817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCO)"""),
)

species(
    label = '[C]#COC(=C)O(1480)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,S} {9,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {6,T}
6 C u1 p0 c0 {5,T}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (359.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,2175,525,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.737649,'amu*angstrom^2'), symmetry=1, barrier=(16.96,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.736369,'amu*angstrom^2'), symmetry=1, barrier=(16.9306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58112,'amu*angstrom^2'), symmetry=1, barrier=(36.3531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30302,0.0624553,-8.9529e-05,5.96994e-08,-1.35322e-11,43289.8,18.7224], Tmin=(100,'K'), Tmax=(715.52,'K')), NASAPolynomial(coeffs=[12.301,0.011403,-4.36955e-06,7.27017e-10,-4.55944e-14,41448.9,-32.5207], Tmin=(715.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
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
    E0 = (89.4539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (526.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (54.7703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (57.6991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (50.7259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (151.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (250.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (295.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (429.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (282.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (553.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (196.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (168.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(O)OC#C(1479)'],
    products = ['HCCO(21)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(85.31,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(5)', '[C]=C(O)OC#C(3930)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(O)OC#C(1479)'],
    products = ['OC1=CC=[C]O1(3931)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0.31, Ea=(50.6264,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH
Exact match found for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(O)OC#C(1479)'],
    products = ['[CH]=C1C=C(O)O1(3932)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['HCCO(21)', 'C#CO(429)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(68.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(7)', 'C#COC#C(3933)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(32.9262,'m^3/(mol*s)'), n=1.67015, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_1CS->C_2CS->C_N-3BrOS->S_N-Sp-2C=1C_Ext-1C-R',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_1CS->C_2CS->C_N-3BrOS->S_N-Sp-2C=1C_Ext-1C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCCO(21)', '[CH]=[C]O(432)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H(22)', '[CH]C(=O)O(3934)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(7)', '[CH]=[C]OC#C(3935)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.54e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', '[CH]C(=O)OC#C(3936)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(9.02696,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[C]#COC(=[CH])O(3937)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.6015e+19,'m^3/(mol*s)'), n=-4.65728, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -23.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)OC#C(1479)'],
    products = ['C#COC([CH2])=O(1013)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#COC(=C)O(1480)'],
    products = ['[CH]=C(O)OC#C(1479)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(544115,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H_TSSD;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1695',
    isomers = [
        '[CH]=C(O)OC#C(1479)',
    ],
    reactants = [
        ('HCCO(21)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1695',
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

