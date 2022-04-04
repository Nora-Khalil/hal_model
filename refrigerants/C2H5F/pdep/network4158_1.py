species(
    label = 'O=[C]OC1(F)[CH]C1(11537)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
6  C u1 p0 c0 {4,S} {5,S} {10,S}
7  C u1 p0 c0 {2,S} {3,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-53.5248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,2950,3150,900,1000,1100,1855,455,950,528.286,528.289,528.29,528.303,528.309,528.309,528.309],'cm^-1')),
        HinderedRotor(inertia=(0.000603982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000603974,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63862,0.0434915,-2.27408e-05,-1.93903e-09,3.47843e-12,-6345.4,20.3882], Tmin=(100,'K'), Tmax=(1108.84,'K')), NASAPolynomial(coeffs=[13.166,0.017199,-7.85843e-06,1.54983e-09,-1.12139e-13,-9841.85,-40.6547], Tmin=(1108.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.5248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdOsH) + ring(Cs-Cs(F)(O2)-Cs) + radical(Csj(Cs-F1sO2sCs)(Cs-CsHH)(H)_ring) + radical((O)CJOC)"""),
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
    label = 'FC1=CC1(5962)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,D} {7,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (101.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,323,467,575,827,1418,1195.53,1197.63,1198.13,2278.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.17,'J/mol'), sigma=(4.94763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=463.62 K, Pc=55.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08997,-0.0083026,0.000111185,-1.92685e-07,1.08455e-10,12179.8,8.10644], Tmin=(10,'K'), Tmax=(559.269,'K')), NASAPolynomial(coeffs=[1.56308,0.0243173,-1.53208e-05,4.6222e-09,-5.34453e-13,12234.9,16.7948], Tmin=(559.269,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(101.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]1[CH]C1(7380)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (273.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,212,367,445,1450,1756.33,1756.66,1757.17,1757.52,1757.63],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2879.02,'J/mol'), sigma=(5.07614,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=449.70 K, Pc=49.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51151,0.00823174,1.42819e-05,-1.69952e-08,5.05698e-12,32904.9,15.0919], Tmin=(100,'K'), Tmax=(1211.99,'K')), NASAPolynomial(coeffs=[4.17684,0.0157606,-7.07148e-06,1.37064e-09,-9.69304e-14,32029.4,8.80755], Tmin=(1211.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsCsF1s) + radical(cyclopropane)"""),
)

species(
    label = '[C]=O(126)',
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
    label = '[O]C1(F)[CH]C1(11606)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {3,S} {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (54.4031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,2950,3150,900,1000,1100,180,540.629,866.039,866.858,867.191,2404.22],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27876,0.0316595,-1.69882e-05,-3.06489e-09,3.95385e-12,6610.59,16.2839], Tmin=(100,'K'), Tmax=(1017.38,'K')), NASAPolynomial(coeffs=[11.1959,0.0095825,-3.57907e-06,6.90842e-10,-5.10726e-14,4124.29,-30.1918], Tmin=(1017.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.4031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCsCs)_1971_ring) + radical(Csj(Cs-F1sO2sCs)(Cs-CsHH)(H)_ring)"""),
)

species(
    label = 'O=C1OC2(F)CC12(11538)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-397.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45132,0.0220278,3.08857e-05,-5.07552e-08,1.87906e-11,-47799.6,13.7262], Tmin=(100,'K'), Tmax=(1052.49,'K')), NASAPolynomial(coeffs=[10.6369,0.0207805,-9.89605e-06,2.03454e-09,-1.51936e-13,-51176.6,-34.0406], Tmin=(1052.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-397.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cds-OdCsOs) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'O=COC1(F)C=C1(9962)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {6,D} {8,S}
6  C u0 p0 c0 {4,S} {5,D} {9,S}
7  C u0 p0 c0 {2,S} {3,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-286.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00638,0.0335791,2.29736e-06,-2.54549e-08,1.10706e-11,-34411.3,19.8413], Tmin=(100,'K'), Tmax=(1062.36,'K')), NASAPolynomial(coeffs=[12.2193,0.0181644,-8.46759e-06,1.71389e-09,-1.26726e-13,-37881.4,-36.167], Tmin=(1062.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-286.782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + ring(Cs(F)(O2)-Cd-Cd)"""),
)

species(
    label = '[CH2]C=C(F)O[C]=O(9073)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,D} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {3,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-167.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,326,540,652,719,1357,3000,3100,440,815,1455,1000,1855,455,950,180,712.385],'cm^-1')),
        HinderedRotor(inertia=(2.04846,'amu*angstrom^2'), symmetry=1, barrier=(47.0981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04712,'amu*angstrom^2'), symmetry=1, barrier=(47.0673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784626,'amu*angstrom^2'), symmetry=1, barrier=(18.0401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57279,0.0575252,-6.59148e-05,4.04516e-08,-1.0234e-11,-20059.6,19.9771], Tmin=(100,'K'), Tmax=(945.734,'K')), NASAPolynomial(coeffs=[9.62178,0.0234814,-1.19184e-05,2.38804e-09,-1.71948e-13,-21582,-18.4055], Tmin=(945.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCFO) + group(Cds-OdOsH) + radical(Allyl_P) + radical((O)CJOC)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.56331e-09,3.13594e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88611e-07,1.21037e-10,-7.84012e-15,-14180.9,6.71043], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'O=[C]OC1=CC1(11607)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {6,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {3,S} {4,D} {9,S}
6 C u1 p0 c0 {1,S} {2,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (128.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,1855,455,950,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73551,0.0427185,-3.18624e-05,6.25938e-09,1.59845e-12,15601,17.4763], Tmin=(100,'K'), Tmax=(1036.41,'K')), NASAPolynomial(coeffs=[13.5178,0.0104609,-4.30358e-06,8.3603e-10,-6.12252e-14,12449,-43.2123], Tmin=(1036.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + ring(Cyclopropene) + radical((O)CJOC)"""),
)

species(
    label = '[O][C]=O(625)',
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
    label = 'O=[C]OC1(F)C=C1(11608)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {7,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u0 p0 c0 {4,S} {5,D} {9,S}
7 C u1 p0 c0 {2,S} {3,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-90.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,1855,455,950,543.285,543.3,543.337,543.382,543.505,543.523],'cm^-1')),
        HinderedRotor(inertia=(0.226578,'amu*angstrom^2'), symmetry=1, barrier=(47.4477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226487,'amu*angstrom^2'), symmetry=1, barrier=(47.4529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76797,0.041659,-2.67802e-05,3.87856e-09,1.35902e-12,-10778,19.3261], Tmin=(100,'K'), Tmax=(1138.84,'K')), NASAPolynomial(coeffs=[12.8644,0.0145782,-6.77688e-06,1.33935e-09,-9.6726e-14,-14076.6,-39.0366], Tmin=(1138.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdOsH) + ring(Cs(F)(O2)-Cd-Cd) + radical((O)CJOC)"""),
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
    label = 'O=[C]OC1=C[CH]1(10004)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {6,D}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 C u1 p0 c0 {1,S} {2,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (328.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,1855,455,950,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0573,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94398,0.0373601,-2.28401e-05,-3.55481e-09,5.48562e-12,39637,18.5304], Tmin=(100,'K'), Tmax=(988.715,'K')), NASAPolynomial(coeffs=[13.8674,0.00652092,-2.44937e-06,4.94472e-10,-3.86085e-14,36428.8,-43.1584], Tmin=(988.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdOsH) + ring(Cyclopropene) + radical(CCJCO) + radical((O)CJOC)"""),
)

species(
    label = 'O=COC1(F)[CH][CH]1(11609)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u1 p0 c0 {4,S} {6,S} {8,S}
6  C u1 p0 c0 {4,S} {5,S} {9,S}
7  C u0 p0 c0 {2,S} {3,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-12.7675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,3150,900,1100,2782.5,750,1395,475,1775,1000,481.835,481.862,481.866,481.89,481.909,481.934],'cm^-1')),
        HinderedRotor(inertia=(0.000725993,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000726022,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72438,0.03933,-7.46575e-06,-1.98396e-08,1.01457e-11,-1444.37,22.0111], Tmin=(100,'K'), Tmax=(1042.83,'K')), NASAPolynomial(coeffs=[13.9526,0.0157814,-7.18792e-06,1.45913e-09,-1.08865e-13,-5264.71,-43.5848], Tmin=(1042.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.7675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdOsH) + ring(Cs-Cs(F)(O2)-Cs) + radical(Csj(Cs-F1sO2sCs)(Cs-CsHH)(H)_ring) + radical(Csj(Cs-F1sO2sCs)(Cs-CsHH)(H)_ring)"""),
)

species(
    label = 'O=[C]O[C]1CC1F(9079)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u1 p0 c0 {2,S} {3,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-49.2087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9282,0.0470313,-4.56819e-05,2.3681e-08,-5.06912e-12,-5845.04,21.7366], Tmin=(100,'K'), Tmax=(1105.9,'K')), NASAPolynomial(coeffs=[9.31676,0.0203076,-9.43557e-06,1.83107e-09,-1.29808e-13,-7479.28,-14.6528], Tmin=(1105.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.2087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cds-OdOsH) + ring(Cs(F)-Cs-Cs(O2)) + radical(C2CsJOC(O)H) + radical((O)CJOCC2)"""),
)

species(
    label = 'O=C(F)O[C]1[CH]C1(11610)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u1 p0 c0 {2,S} {4,S} {6,S}
6  C u1 p0 c0 {4,S} {5,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-103.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93459,0.0357369,-6.94568e-06,-1.83221e-08,9.55165e-12,-12417,25.0769], Tmin=(100,'K'), Tmax=(1029.1,'K')), NASAPolynomial(coeffs=[13.099,0.0136984,-5.95171e-06,1.19983e-09,-8.9765e-14,-15845.8,-34.5995], Tmin=(1029.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(COFOO) + ring(Cyclopropane) + radical(C2CsJOC(O)) + radical(CCJCO)"""),
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
    E0 = (-61.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (485.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-53.4737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-36.7848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (6.07656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-31.8661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-5.88931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (237.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (129.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (113.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (296.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (254.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (71.2762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (128.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (140.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]OC1(F)[CH]C1(11537)'],
    products = ['CO2(14)', 'FC1=CC1(5962)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=O(126)', '[O]C1(F)[CH]C1(11606)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]OC1(F)[CH]C1(11537)'],
    products = ['O=C1OC2(F)CC12(11538)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]OC1(F)[CH]C1(11537)'],
    products = ['O=COC1(F)C=C1(9962)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C=C(F)O[C]=O(9073)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3_D;doublebond_intra_pri;radadd_intra_cs] + [R3_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CO(13)', '[O]C1(F)[CH]C1(11606)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(34.1,'m^3/(mol*s)'), n=8.73864e-09, Ea=(40.7048,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CO2(14)', 'F[C]1[CH]C1(7380)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.33888e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(132.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'O=[C]OC1=CC1(11607)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(43.6851,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=O(625)', 'FC1=CC1(5962)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(4.61629,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'O=[C]OC1(F)C=C1(11608)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.07428,'m^3/(mol*s)'), n=2.32905, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.015350652311676958, var=0.006503997552278654, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(625)', 'F[C]1[CH]C1(7380)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['HF(38)', 'O=[C]OC1=C[CH]1(10004)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(562.232,'m^3/(mol*s)'), n=1.03051, Ea=(214.591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=COC1(F)[CH][CH]1(11609)'],
    products = ['O=[C]OC1(F)[CH]C1(11537)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9118.87,'s^-1'), n=2.33789, Ea=(92.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;XH_out] for rate rule [R4H_SSS;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]OC1(F)[CH]C1(11537)'],
    products = ['O=[C]O[C]1CC1F(9079)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(190.747,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]OC1(F)[CH]C1(11537)'],
    products = ['O=C(F)O[C]1[CH]C1(11610)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(202.628,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #4158',
    isomers = [
        'O=[C]OC1(F)[CH]C1(11537)',
    ],
    reactants = [
        ('CO2(14)', 'FC1=CC1(5962)'),
        ('CO2(14)', 'F[C]1[CH]C1(7380)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4158',
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

