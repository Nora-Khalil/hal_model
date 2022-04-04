species(
    label = '[O]C(=O)O[CH]O(754)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {2,S} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-366.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4292.31,'J/mol'), sigma=(6.39584,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=670.45 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75959,0.0454949,-5.16469e-05,2.81572e-08,-5.90453e-12,-44020.1,20.8387], Tmin=(100,'K'), Tmax=(1177.28,'K')), NASAPolynomial(coeffs=[13.0314,0.00719721,-2.85097e-06,5.25275e-10,-3.67809e-14,-46674.1,-35.3808], Tmin=(1177.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(OCJO)"""),
)

species(
    label = 'O=COC(=O)O(758)',
    structure = adjacencyList("""1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 C u0 p0 c0 {1,S} {4,D} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-722.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,456.563,456.583,456.745,456.78,457.058,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000807833,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145867,'amu*angstrom^2'), symmetry=1, barrier=(21.6117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14594,'amu*angstrom^2'), symmetry=1, barrier=(21.6078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4118.47,'J/mol'), sigma=(5.92789,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=643.30 K, Pc=44.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01007,0.0367916,-2.7478e-05,4.26228e-09,1.91234e-12,-86767.4,21.1977], Tmin=(100,'K'), Tmax=(1050.58,'K')), NASAPolynomial(coeffs=[13.3571,0.00666628,-3.13751e-06,6.65153e-10,-5.12333e-14,-89873.3,-37.5397], Tmin=(1050.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-722.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs) + group(Cds-OdOsH)"""),
)

species(
    label = '[O]C(O)O[C]=O(755)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u1 p0 c0 {1,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-327.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,1855,455,950,180,1048.97,1049.13],'cm^-1')),
        HinderedRotor(inertia=(0.00366082,'amu*angstrom^2'), symmetry=1, barrier=(2.86089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124115,'amu*angstrom^2'), symmetry=1, barrier=(2.85365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124369,'amu*angstrom^2'), symmetry=1, barrier=(2.85949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4285.1,'J/mol'), sigma=(6.64969,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.32 K, Pc=33.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5317,0.0485168,-0.000111882,1.35217e-07,-5.71213e-11,-39400.7,19.7659], Tmin=(100,'K'), Tmax=(839.773,'K')), NASAPolynomial(coeffs=[-4.87974,0.0413573,-2.32487e-05,4.64392e-09,-3.25143e-13,-36658.7,63.1421], Tmin=(839.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH) + radical(OCOJ) + radical((O)CJOC)"""),
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
    label = 'O=CO(77)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-389.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,910.228,910.24],'cm^-1')),
        HinderedRotor(inertia=(0.000203447,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76624,-0.000867782,3.32891e-05,-3.95828e-08,1.39946e-11,-46830.2,7.56886], Tmin=(100,'K'), Tmax=(997.616,'K')), NASAPolynomial(coeffs=[5.52461,0.0077913,-3.35088e-06,6.86825e-10,-5.23916e-14,-47962.8,-4.82803], Tmin=(997.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'O=C(O)O(747)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,421.686,423.663,424.543,424.784],'cm^-1')),
        HinderedRotor(inertia=(0.00093875,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344145,'amu*angstrom^2'), symmetry=1, barrier=(43.4791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.87,'J/mol'), sigma=(5.42485,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.51 K, Pc=51.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98124,0.0114379,2.77981e-05,-4.94746e-08,2.0802e-11,-75136.3,10.9454], Tmin=(100,'K'), Tmax=(956.736,'K')), NASAPolynomial(coeffs=[11.8951,0.00171649,-1.48201e-07,9.26559e-11,-1.38667e-14,-78102.6,-38.253], Tmin=(956.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs)"""),
)

species(
    label = 'O=[C]O(125)',
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
    label = '[O]C=O(121)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-138.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51545,0.0112656,-9.5808e-06,4.36315e-09,-8.44574e-13,-16672.2,7.37034], Tmin=(100,'K'), Tmax=(1184.09,'K')), NASAPolynomial(coeffs=[5.09942,0.00591471,-2.80223e-06,5.46638e-10,-3.87735e-14,-17047.4,-0.539013], Tmin=(1184.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""formyloxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HCO(16)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C(=O)O(309)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-381.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,625.974,626.005,626.803,627.197,627.93],'cm^-1')),
        HinderedRotor(inertia=(0.000430732,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.87,'J/mol'), sigma=(5.42485,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.51 K, Pc=51.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28603,0.00907564,1.53519e-05,-2.86387e-08,1.1972e-11,-45835.1,11.0209], Tmin=(100,'K'), Tmax=(968.183,'K')), NASAPolynomial(coeffs=[8.68687,0.00322549,-1.09079e-06,2.46308e-10,-2.1583e-14,-47652.5,-18.8452], Tmin=(968.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-381.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs) + radical(OC=OOJ)"""),
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
    label = 'O=[C]OC=O(230)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {6,S}
5 C u1 p0 c0 {1,S} {3,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-288.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1855,455,950,378.563],'cm^-1')),
        HinderedRotor(inertia=(0.0786445,'amu*angstrom^2'), symmetry=1, barrier=(7.99249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0785762,'amu*angstrom^2'), symmetry=1, barrier=(7.99332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0274,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.87,'J/mol'), sigma=(5.46152,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.17 K, Pc=48.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8329,0.027821,-3.45944e-05,2.30648e-08,-6.29073e-12,-34627,15.953], Tmin=(100,'K'), Tmax=(884.048,'K')), NASAPolynomial(coeffs=[6.62732,0.0106531,-5.46551e-06,1.09902e-09,-7.91882e-14,-35297.9,-1.88535], Tmin=(884.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-OdOsH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[O]C([O])O(220)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {6,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (-145.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2358.11,2358.29],'cm^-1')),
        HinderedRotor(inertia=(3.03187e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4214.4,'J/mol'), sigma=(6.62052,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.28 K, Pc=32.95 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79685,0.0288372,-0.000108613,1.59396e-07,-7.20139e-11,-17564.9,12.6702], Tmin=(100,'K'), Tmax=(855.716,'K')), NASAPolynomial(coeffs=[-15.1237,0.0502966,-2.88126e-05,5.74781e-09,-3.9998e-13,-11874.4,115.333], Tmin=(855.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[O]C([O])=O(268)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-60.3782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,620.879,889.717,889.721,889.771,1871.84],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.87,'J/mol'), sigma=(5.42485,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.51 K, Pc=51.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36989,0.0105617,-1.2751e-06,-7.37812e-09,3.82317e-12,-7236.21,10.2409], Tmin=(100,'K'), Tmax=(995.957,'K')), NASAPolynomial(coeffs=[7.17797,0.00288732,-1.19276e-06,2.48511e-10,-1.9462e-14,-8372.66,-10.0127], Tmin=(995.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.3782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CO3t1""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]O(118)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.873,3308.85],'cm^-1')),
        HinderedRotor(inertia=(0.0103145,'amu*angstrom^2'), symmetry=1, barrier=(22.7124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75937,0.00296141,8.9037e-06,-1.3501e-08,5.39791e-12,24775.6,6.76289], Tmin=(100,'K'), Tmax=(940.438,'K')), NASAPolynomial(coeffs=[5.09115,0.00321234,-9.31657e-07,1.59608e-10,-1.15723e-14,24263.5,-0.97116], Tmin=(940.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'O=[C]O[CH]O(274)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {6,S}
5 C u1 p0 c0 {1,S} {3,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-176.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1855,455,950,446.436,448.147],'cm^-1')),
        HinderedRotor(inertia=(0.209429,'amu*angstrom^2'), symmetry=1, barrier=(29.6476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0648234,'amu*angstrom^2'), symmetry=1, barrier=(9.17277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.560735,'amu*angstrom^2'), symmetry=1, barrier=(79.832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0354,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17841,0.0401061,-5.0797e-05,3.21969e-08,-7.96432e-12,-21176.7,16.5643], Tmin=(100,'K'), Tmax=(994.594,'K')), NASAPolynomial(coeffs=[9.92545,0.00895057,-3.81111e-06,7.03717e-10,-4.84946e-14,-22717.7,-20.7689], Tmin=(994.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-OdOsH) + radical(OCJO) + radical((O)CJOC)"""),
)

species(
    label = 'O=C1OC(O)O1(762)',
    structure = adjacencyList("""1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-643.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10159,0.00571088,5.76784e-05,-7.59603e-08,2.76843e-11,-77346.1,16.3978], Tmin=(100,'K'), Tmax=(1012.55,'K')), NASAPolynomial(coeffs=[10.4431,0.0128136,-6.33021e-06,1.39932e-09,-1.109e-13,-80683.7,-28.2519], Tmin=(1012.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsOs) + ring(Cyclobutane)"""),
)

species(
    label = 'O[CH]O[C]1OO1(766)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {3,S}
6 C u1 p0 c0 {1,S} {4,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (16.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7362,0.0495162,-6.74959e-05,4.57852e-08,-1.19752e-11,2066.6,21.0931], Tmin=(100,'K'), Tmax=(947.717,'K')), NASAPolynomial(coeffs=[11.3718,0.00884839,-3.1305e-06,5.08854e-10,-3.19791e-14,240.184,-24.8761], Tmin=(947.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + ring(dioxirane) + radical(Cs_P) + radical(OCJO)"""),
)

species(
    label = '[O][C]1OC(O)O1(767)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u1 p0 c0 {1,S} {2,S} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-180.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.63983,0.00460057,1.85914e-05,-1.51027e-08,2.92848e-12,-21772.5,9.49065], Tmin=(100,'K'), Tmax=(2017.39,'K')), NASAPolynomial(coeffs=[15.7323,0.014902,-1.16063e-05,2.32437e-09,-1.54089e-13,-32013,-61.5707], Tmin=(2017.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsOsH) + ring(Cyclobutane) + radical(OCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C1([O])OC1O(768)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-196.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38681,0.0425718,-7.15334e-05,6.64077e-08,-2.23881e-11,-23574.3,19.8524], Tmin=(100,'K'), Tmax=(968.09,'K')), NASAPolynomial(coeffs=[2.77021,0.0211312,-7.54603e-06,1.15638e-09,-6.62292e-14,-22718.1,22.8209], Tmin=(968.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-CsOsOsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]C(=O)OC=O(616)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-478.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,180,180,897.583,897.642,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0589541,'amu*angstrom^2'), symmetry=1, barrier=(18.2187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318629,'amu*angstrom^2'), symmetry=1, barrier=(18.2189,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0268,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47932,0.0325023,-3.32199e-05,1.64079e-08,-3.195e-12,-57473.4,19.9889], Tmin=(100,'K'), Tmax=(1239.84,'K')), NASAPolynomial(coeffs=[9.87762,0.0086336,-4.34259e-06,8.80426e-10,-6.40333e-14,-59308,-17.2941], Tmin=(1239.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-478.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs) + group(Cds-OdOsH) + radical(OC=OOJ)"""),
)

species(
    label = '[O][CH]O(212)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1989.6],'cm^-1')),
        HinderedRotor(inertia=(0.0936924,'amu*angstrom^2'), symmetry=1, barrier=(2.15417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3630.86,'J/mol'), sigma=(5.89779,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.13 K, Pc=40.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40205,0.0209569,-4.96195e-05,5.95345e-08,-2.44629e-11,660.882,10.3083], Tmin=(100,'K'), Tmax=(868.742,'K')), NASAPolynomial(coeffs=[-0.343261,0.0179326,-9.40006e-06,1.81362e-09,-1.23836e-13,2076.49,32.2525], Tmin=(868.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = '[O]COC([O])=O(769)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {3,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-328.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.36424,0.0262322,4.98969e-06,-2.75967e-08,1.25315e-11,-39413.5,19.3551], Tmin=(100,'K'), Tmax=(1008.88,'K')), NASAPolynomial(coeffs=[12.2605,0.00974511,-4.32114e-06,9.06774e-10,-7.02642e-14,-42568.1,-34.2138], Tmin=(1008.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(OCOJ) + radical(OC=OOJ)"""),
)

species(
    label = '[O][CH]OC(=O)O(764)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 C u1 p0 c0 {1,S} {4,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-383.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39903,0.0400992,-4.49791e-05,3.03493e-08,-9.18082e-12,-46017.2,18.5626], Tmin=(100,'K'), Tmax=(769.349,'K')), NASAPolynomial(coeffs=[5.54618,0.0237367,-1.30774e-05,2.70567e-09,-1.98102e-13,-46501.4,4.20456], Tmin=(769.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = 'O[C]1O[CH]OO1(810)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {5,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u1 p0 c0 {1,S} {3,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-46.6465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654622,0.0444707,-4.51207e-05,2.19331e-08,-3.74843e-12,-5464.22,24.9902], Tmin=(100,'K'), Tmax=(1857.49,'K')), NASAPolynomial(coeffs=[9.86668,0.00524434,2.21333e-06,-6.74887e-10,5.06914e-14,-5541.66,-16.1534], Tmin=(1857.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.6465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + ring(124trioxolane) + radical(Cs_P) + radical(OCJO)"""),
)

species(
    label = '[O]C([O])OC=O(237)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-296.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1233.12,1233.15,1233.17,1233.2],'cm^-1')),
        HinderedRotor(inertia=(0.00329141,'amu*angstrom^2'), symmetry=1, barrier=(3.55194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00329166,'amu*angstrom^2'), symmetry=1, barrier=(3.55203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[6.92093,-0.00339225,3.1511e-05,-2.18273e-08,4.00675e-12,-35875.9,5.04701], Tmin=(100,'K'), Tmax=(2057.32,'K')), NASAPolynomial(coeffs=[17.8531,0.0171257,-1.39057e-05,2.7593e-09,-1.80262e-13,-49214.4,-77.0658], Tmin=(2057.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-296.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = 'O=[C]O[C](O)O(763)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {2,S} {3,S}
6 C u1 p0 c0 {1,S} {4,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-350.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,1855,455,950,258.851,258.946],'cm^-1')),
        HinderedRotor(inertia=(0.00251619,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00251463,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00251627,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00251491,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0348,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01929,0.0519935,-9.72831e-05,9.69217e-08,-3.6558e-11,-42038.8,21.129], Tmin=(100,'K'), Tmax=(833.414,'K')), NASAPolynomial(coeffs=[3.98522,0.0238583,-1.29883e-05,2.57008e-09,-1.79335e-13,-41717.1,15.8989], Tmin=(833.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-350.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH) + radical(Cs_P) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OC(=O)O(632)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 C u1 p0 c0 {1,S} {4,D}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-525.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,491.58,491.599,491.599,491.6,491.608,491.619],'cm^-1')),
        HinderedRotor(inertia=(0.0153872,'amu*angstrom^2'), symmetry=1, barrier=(2.6384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0793718,'amu*angstrom^2'), symmetry=1, barrier=(13.6143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0793872,'amu*angstrom^2'), symmetry=1, barrier=(13.6142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0268,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4118.47,'J/mol'), sigma=(5.92789,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=643.30 K, Pc=44.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85324,0.0439492,-5.35247e-05,2.99567e-08,-6.37671e-12,-63137.7,20.3873], Tmin=(100,'K'), Tmax=(1164.52,'K')), NASAPolynomial(coeffs=[13.526,0.00385427,-1.87881e-06,3.90192e-10,-2.93403e-14,-65856.3,-37.7051], Tmin=(1164.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-525.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    E0 = (-169.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (342.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (263.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-160.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-144.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (214.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (17.2072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1.00107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-95.1947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-40.0336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-127.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (234.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-27.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-75.6869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-292.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (207.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-76.6554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-67.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-118.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-105.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-133.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-151.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-62.4067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-69.0799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-116.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-403.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-130.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (263.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-122.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-54.7061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-122.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-87.3023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-19.1771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-71.7292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (234.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (27.5702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (34.3157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (27.8768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['CO2(15)', 'O=CO(77)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C([O])=O(268)', '[CH]O(118)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(131316,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(7)', 'O=[C]O[CH]O(274)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['O=C1OC(O)O1(762)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['O[CH]O[C]1OO1(766)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(383.716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['[O][C]1OC(O)O1(767)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.29795e+08,'s^-1'), n=0.935645, Ea=(186.471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 1.7320508075688772
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 184.0 to 186.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['[O]C1([O])OC1O(768)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.44049e+09,'s^-1'), n=0.679905, Ea=(170.265,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHNd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 169.6 to 170.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=O(213)', 'O=CO(77)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(65.3144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', '[O]C(=O)OC=O(616)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.07,'m^3/(mol*s)'), n=2.43, Ea=(29.0463,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO2(15)', '[O][CH]O(212)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(73.0345,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(213)', '[O][CH]O(212)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]COC([O])=O(769)'],
    products = ['[O]C(=O)O[CH]O(754)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.41377e+07,'s^-1'), n=1.57383, Ea=(103.529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(=O)O[CH]O(754)'],
    products = ['[O][CH]OC(=O)O(764)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(72901.1,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_3;O_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(14)', 'O=C(O)O(747)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.54e-07,'m^3/(mol*s)'), n=3.7, Ea=(254.351,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs_2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->H_1COCbCdCtHNOSSidSis->O',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs_2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->H_1COCbCdCtHNOSSidSis->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O[C]1O[CH]OO1(810)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.43612e+09,'s^-1'), n=0.0758668, Ea=(56.4791,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C([O])OC=O(237)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(O)O[C]=O(755)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]O[C](O)O(763)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]COC([O])=O(769)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]O(125)', '[O]C=O(121)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['HCO(16)', '[O]C(=O)O(309)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OH(8)', 'O=[C]OC=O(230)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(6)', '[O]C(=O)OC=O(616)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(6)', 'O=[C]OC(=O)O(632)'],
    products = ['O=COC(=O)O(758)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.31504e+12,'m^3/(mol*s)'), n=-2.20453, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O"""),
)

reaction(
    label = 'reaction1',
    reactants = ['O=COC(=O)O(758)'],
    products = ['CO2(15)', 'O=CO(77)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.82035e+06,'s^-1'), n=1.67734, Ea=(120.903,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2285837836968013, var=2.96975919141528, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(O)O[C]=O(755)'],
    products = ['CO2(15)', 'O=CO(77)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(7)', 'O=[C]O[CH]O(274)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/O2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(O)O[C]=O(755)'],
    products = ['O=C1OC(O)O1(762)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CO(14)', '[O]C([O])O(220)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(68.2,'m^3/(mol*s)'), n=8.73864e-09, Ea=(12.5045,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CO2(15)', '[O][CH]O(212)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(122.376,'m^3/(mol*s)'), n=1.29312, Ea=(77.8224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.28976280499384915, var=2.1569028208455543, Tref=1000.0, N=37, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]=O(213)', 'O=CO(77)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(73.2068,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(8)', 'O=[C]OC=O(230)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(43.2296,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', 'O=[C]OC(=O)O(632)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(44.6505,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(213)', '[O][CH]O(212)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(O)O[C]=O(755)'],
    products = ['O=[C]O[C](O)O(763)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(O)O[C]=O(755)'],
    products = ['[O][CH]OC(=O)O(764)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.9172e+08,'s^-1'), n=1.32036, Ea=(164.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C([O])OC=O(237)'],
    products = ['[O]C(O)O[C]=O(755)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.50344e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #404',
    isomers = [
        '[O]C(=O)O[CH]O(754)',
        'O=COC(=O)O(758)',
        '[O]C(O)O[C]=O(755)',
    ],
    reactants = [
        ('CO2(15)', 'O=CO(77)'),
        ('[O][C]=O(213)', 'O=CO(77)'),
        ('CO(14)', 'O=C(O)O(747)'),
        ('O=[C]O(125)', '[O]C=O(121)'),
        ('HCO(16)', '[O]C(=O)O(309)'),
        ('OH(8)', 'O=[C]OC=O(230)'),
        ('CO(14)', '[O]C([O])O(220)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #404',
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

