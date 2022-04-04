species(
    label = '[CH2]CC=[C]CF(11807)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (215.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,122,584,932,1043,1247,1392,1461,2986,3039,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,301.621],'cm^-1')),
        HinderedRotor(inertia=(0.0946286,'amu*angstrom^2'), symmetry=1, barrier=(6.11278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947148,'amu*angstrom^2'), symmetry=1, barrier=(6.11278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0134113,'amu*angstrom^2'), symmetry=1, barrier=(21.2824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75355,0.0452342,-2.69131e-05,7.59206e-09,-8.53632e-13,25983.8,25.3257], Tmin=(100,'K'), Tmax=(2013.23,'K')), NASAPolynomial(coeffs=[14.1563,0.0205918,-8.55302e-06,1.51229e-09,-9.86625e-14,20989.8,-43.1894], Tmin=(2013.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = 'C2H4(30)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C#CCF(3587)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (8.12032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([319,1023,1071,1259,1317,1409,3054,3019,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.60388,'amu*angstrom^2'), symmetry=1, barrier=(36.8763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2742,'J/mol'), sigma=(4.732,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.29 K, Pc=58.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95512,0.00280826,7.03906e-05,-1.43186e-07,9.06266e-11,979.285,8.14914], Tmin=(10,'K'), Tmax=(503.334,'K')), NASAPolynomial(coeffs=[2.71827,0.0217978,-1.34993e-05,4.0832e-09,-4.79139e-13,987.76,12.1145], Tmin=(503.334,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(8.12032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""C#CCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C=[C]CF(11864)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {5,D} {8,S}
4  C u1 p0 c0 {3,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {3,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (182.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.177498,'amu*angstrom^2'), symmetry=1, barrier=(4.08103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93869,'amu*angstrom^2'), symmetry=1, barrier=(44.5743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55728,0.0328081,-1.90353e-05,5.25907e-09,-5.85398e-13,22001.2,17.4977], Tmin=(100,'K'), Tmax=(1959.86,'K')), NASAPolynomial(coeffs=[9.80254,0.018021,-7.71788e-06,1.40938e-09,-9.43358e-14,19161.2,-22.3317], Tmin=(1959.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Allyl_P) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = 'FCC1=CCC1(12011)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-70.0879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14474,0.0276093,3.30749e-05,-5.87113e-08,2.30553e-11,-8351.18,18.3662], Tmin=(100,'K'), Tmax=(999.676,'K')), NASAPolynomial(coeffs=[11.1595,0.023721,-9.38018e-06,1.80469e-09,-1.32099e-13,-11761.6,-33.1651], Tmin=(999.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.0879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=CC=CCF(12013)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-117.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57663,0.0441518,-1.25985e-05,-1.40607e-08,8.09515e-12,-14049.4,20.1495], Tmin=(100,'K'), Tmax=(1031.26,'K')), NASAPolynomial(coeffs=[12.2194,0.0221973,-8.7755e-06,1.64022e-09,-1.16448e-13,-17272.2,-36.5058], Tmin=(1031.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CCC=C=CF(8207)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {1,S} {6,D} {13,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-49.7003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.606,0.0459216,-2.60873e-05,5.5609e-09,4.41051e-14,-5885.86,21.0515], Tmin=(100,'K'), Tmax=(1286.93,'K')), NASAPolynomial(coeffs=[10.9953,0.0243884,-9.90569e-06,1.79758e-09,-1.22167e-13,-8936.03,-29.0762], Tmin=(1286.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.7003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=[C]CF(4166)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u1 p0 c0 {2,S} {4,D}
4 C u1 p0 c0 {3,D} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (311.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,380.029,380.029,380.029,1901.95,1901.95,1901.95,1901.95,3346.47],'cm^-1')),
        HinderedRotor(inertia=(0.117101,'amu*angstrom^2'), symmetry=1, barrier=(12.0011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74483,0.0317272,-5.02121e-05,4.7676e-08,-1.76456e-11,37561.4,13.6835], Tmin=(100,'K'), Tmax=(834.53,'K')), NASAPolynomial(coeffs=[4.04639,0.0168537,-7.95765e-06,1.52225e-09,-1.05087e-13,37644.8,9.44123], Tmin=(834.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_S) + radical(Cds_P)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,4.63019e-14,-6.5121e-17,3.00122e-20,-4.26132e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3821.96,'K')), NASAPolynomial(coeffs=[2.5,2.03348e-10,-7.42469e-14,1.19914e-17,-7.22693e-22,25474.2,-0.444972], Tmin=(3821.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=CC=[C]CF(11879)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {5,D} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u1 p0 c0 {2,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (120.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.829807,'amu*angstrom^2'), symmetry=1, barrier=(19.0789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828553,'amu*angstrom^2'), symmetry=1, barrier=(19.0501,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51619,0.0489796,-3.76904e-05,1.48269e-08,-2.35482e-12,14555.2,20.8459], Tmin=(100,'K'), Tmax=(1484.81,'K')), NASAPolynomial(coeffs=[12.5005,0.0193886,-7.7967e-06,1.40489e-09,-9.49475e-14,11293.3,-36.489], Tmin=(1484.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = '[CH2]CC=C=C(6661)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {11,S} {12,S}
5  C u0 p0 c0 {2,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (329.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,597.038],'cm^-1')),
        HinderedRotor(inertia=(0.41105,'amu*angstrom^2'), symmetry=1, barrier=(9.45084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818849,'amu*angstrom^2'), symmetry=1, barrier=(18.8269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02941,0.0367789,-1.08451e-05,-8.05917e-09,4.61377e-12,39682.5,19.8002], Tmin=(100,'K'), Tmax=(1084.8,'K')), NASAPolynomial(coeffs=[9.2456,0.022862,-9.1506e-06,1.68413e-09,-1.17036e-13,37370.1,-19.043], Tmin=(1084.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=C=CF(8208)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {6,D} {12,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (155.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,113,247,382,1207,3490,540,610,2055,860.323,864.443],'cm^-1')),
        HinderedRotor(inertia=(0.00731703,'amu*angstrom^2'), symmetry=1, barrier=(3.88244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166865,'amu*angstrom^2'), symmetry=1, barrier=(3.83655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80549,0.0448486,-3.0645e-05,1.06659e-08,-1.52529e-12,18789.4,22.114], Tmin=(100,'K'), Tmax=(1595.44,'K')), NASAPolynomial(coeffs=[10.9744,0.0218607,-9.03232e-06,1.63481e-09,-1.10152e-13,15863.8,-26.4038], Tmin=(1595.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH2](611)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,971.403,1437.07,2688.5,2689.57],'cm^-1')),
        HinderedRotor(inertia=(0.0117751,'amu*angstrom^2'), symmetry=1, barrier=(17.243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40737,0.0100311,6.40964e-06,-1.41296e-08,5.92694e-12,38288.2,6.11701], Tmin=(100,'K'), Tmax=(954.253,'K')), NASAPolynomial(coeffs=[5.52246,0.00856177,-2.90746e-06,5.0236e-10,-3.44577e-14,37547.8,-5.75261], Tmin=(954.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]CC#CCF(12014)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {6,T}
6  C u0 p0 c0 {3,S} {5,T}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (137.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,319,1023,1071,1259,1317,1409,3054,3019,3000,3100,440,815,1455,1000,2100,2250,500,550,299.55,299.602],'cm^-1')),
        HinderedRotor(inertia=(0.0833963,'amu*angstrom^2'), symmetry=1, barrier=(5.31134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08499,'amu*angstrom^2'), symmetry=1, barrier=(69.127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256424,'amu*angstrom^2'), symmetry=1, barrier=(69.1266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78235,0.0443865,-3.24584e-05,1.3684e-08,-2.42932e-12,16578.4,22.7369], Tmin=(100,'K'), Tmax=(1317.83,'K')), NASAPolynomial(coeffs=[8.85257,0.0229262,-8.03152e-06,1.32684e-09,-8.50849e-14,14714.9,-13.3241], Tmin=(1317.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C=[C]CF(12015)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {5,S} {12,S}
5  C u0 p0 c0 {4,S} {6,D} {13,S}
6  C u1 p0 c0 {3,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (151.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,122,584,932,1043,1247,1392,1461,2986,3039,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.0107321,'amu*angstrom^2'), symmetry=1, barrier=(8.23482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04676,'amu*angstrom^2'), symmetry=1, barrier=(24.067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106033,'amu*angstrom^2'), symmetry=1, barrier=(24.0555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84893,0.0401749,-1.12211e-05,-7.58306e-09,3.90364e-12,18269.4,22.0816], Tmin=(100,'K'), Tmax=(1190.88,'K')), NASAPolynomial(coeffs=[9.8223,0.0269393,-1.16117e-05,2.18703e-09,-1.52504e-13,15409.8,-21.8111], Tmin=(1190.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = '[CH2]CC=C[CH]F(8212)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {6,S} {10,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {4,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (108.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,234,589,736,816,1240,3237,953.94],'cm^-1')),
        HinderedRotor(inertia=(0.184563,'amu*angstrom^2'), symmetry=1, barrier=(4.24348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184853,'amu*angstrom^2'), symmetry=1, barrier=(4.25013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59443,'amu*angstrom^2'), symmetry=1, barrier=(36.659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43311,0.0488627,-3.16365e-05,9.97279e-09,-1.25731e-12,13143.3,24.0305], Tmin=(100,'K'), Tmax=(1835.24,'K')), NASAPolynomial(coeffs=[14.6292,0.0201008,-8.12821e-06,1.43313e-09,-9.40081e-14,8299.75,-47.6453], Tmin=(1835.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Csj(Cd-CdH)(F1s)(H))"""),
)

species(
    label = '[CH2]C[C]=CCF(12016)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  C u1 p0 c0 {2,S} {4,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (215.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,212,931,1104,1251,1325,1474,3102,3155,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,352.732,352.767],'cm^-1')),
        HinderedRotor(inertia=(0.00135479,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0826762,'amu*angstrom^2'), symmetry=1, barrier=(7.29901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0826132,'amu*angstrom^2'), symmetry=1, barrier=(7.29752,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09944,0.0463691,-3.1814e-05,1.20121e-08,-2.06079e-12,25927.8,22.6147], Tmin=(100,'K'), Tmax=(1228.42,'K')), NASAPolynomial(coeffs=[6.33185,0.0325876,-1.49858e-05,2.87947e-09,-2.02187e-13,24887.9,1.32499], Tmin=(1228.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'CC[C]=[C]CF(12017)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {6,D}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
"""),
    E0 = (247.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,122,584,932,1043,1247,1392,1461,2986,3039,1670,1700,300,440,280.645],'cm^-1')),
        HinderedRotor(inertia=(0.00222643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192427,'amu*angstrom^2'), symmetry=1, barrier=(10.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185109,'amu*angstrom^2'), symmetry=1, barrier=(10.5344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13502,0.0442821,-2.64005e-05,7.62323e-09,-9.14227e-13,29883.2,22.6052], Tmin=(100,'K'), Tmax=(1757.75,'K')), NASAPolynomial(coeffs=[9.37863,0.0277982,-1.23336e-05,2.28801e-09,-1.55408e-13,27336.8,-16.4266], Tmin=(1757.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = '[CH2]C=C[CH]CF(12018)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u1 p0 c0 {2,S} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u1 p0 c0 {5,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (51.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65431,0.0420283,-5.987e-06,-1.94409e-08,9.4728e-12,6254.1,20.4183], Tmin=(100,'K'), Tmax=(1044.29,'K')), NASAPolynomial(coeffs=[11.6593,0.0242256,-9.88993e-06,1.86698e-09,-1.32801e-13,3045.6,-33.6406], Tmin=(1044.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'CCC=[C][CH]F(8214)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u1 p0 c0 {1,S} {6,S} {13,S}
6  C u1 p0 c0 {4,D} {5,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (141.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6132,0.0451071,-2.07841e-05,-6.19219e-10,2.11611e-12,17093,23.5038], Tmin=(100,'K'), Tmax=(1200.83,'K')), NASAPolynomial(coeffs=[11.1297,0.0252386,-1.07448e-05,2.01215e-09,-1.399e-13,13954.4,-27.7014], Tmin=(1200.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = '[CH2][C]=CCCF(6422)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u1 p0 c0 {6,S} {12,S} {13,S}
6  C u1 p0 c0 {4,D} {5,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (147.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4962,0.0509423,-3.62185e-05,1.31501e-08,-1.95611e-12,17888.5,22.5682], Tmin=(100,'K'), Tmax=(1547.84,'K')), NASAPolynomial(coeffs=[12.0295,0.0237217,-9.83917e-06,1.78826e-09,-1.20996e-13,14627.8,-32.8502], Tmin=(1547.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC=C([CH2])F(8231)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {1,S} {3,D} {6,S}
5  C u1 p0 c0 {2,S} {10,S} {11,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (90.1325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36602,0.0523821,-3.84623e-05,1.45e-08,-2.22562e-12,10939.7,23.7996], Tmin=(100,'K'), Tmax=(1519.47,'K')), NASAPolynomial(coeffs=[12.6384,0.0227073,-9.16739e-06,1.64677e-09,-1.10835e-13,7514.15,-35.2988], Tmin=(1519.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.1325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCsCdF) + radical(RCCJ) + radical(Csj(Cd-F1sCd)(H)(H))"""),
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
    E0 = (36.3184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (384.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (44.6027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (99.7186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (61.2917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (174.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (153.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (285.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (192.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (252.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (195.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (451.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (172.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (175.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (270.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (220.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (186.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (91.3691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (202.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (135.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['C2H4(30)', 'C#CCF(3587)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(18)', '[CH2]C=[C]CF(11864)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['FCC1=CCC1(12011)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['C=CC=CCF(12013)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['CCC=C=CF(8207)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H4(30)', '[CH]=[C]CF(4166)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(108.514,'m^3/(mol*s)'), n=1.59012, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_N-Sp-4R!H-3R_Ext-5R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_N-Sp-4R!H-3R_Ext-5R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'C=CC=[C]CF(11879)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000216213,'m^3/(mol*s)'), n=3.29826, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.02158166639280058, var=0.8466995369935014, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_N-Sp-5R!H-4CClOS_Sp-4CClOS-1CCClOOSS_Sp-2CS=1CCOSS_Sp-5R!H=4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_N-Sp-5R!H-4CClOS_Sp-4CClOS-1CCClOOSS_Sp-2CS=1CCOSS_Sp-5R!H=4CClOS"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', '[CH2]CC=C=C(6661)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.5387,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH2]CC=C=CF(8208)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.296016,'m^3/(mol*s)'), n=2.36777, Ea=(4.20178,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1084748591603159, var=0.5872085620069725, Tref=1000.0, N=10, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH2](611)', 'C#CCF(3587)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.75673e+24,'m^3/(mol*s)'), n=-5.41269, Ea=(105.316,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=2.33529332423437, var=201.40057549922756, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R_4R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[CH2]CC#CCF(12014)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6900,'m^3/(mol*s)'), n=1.64, Ea=(25.2341,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_N-Sp-7R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_N-Sp-7R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH2](611)', '[CH]=[C]CF(4166)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.63958e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]C=[C]CF(12015)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.18e+09,'s^-1'), n=1.22, Ea=(199.995,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 157 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['[CH2]CC=C[CH]F(8212)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C[C]=CCF(12016)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC[C]=[C]CF(12017)'],
    products = ['[CH2]CC=[C]CF(11807)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['[CH2]C=C[CH]CF(12018)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['CCC=[C][CH]F(8214)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(126713,'s^-1'), n=1.75034, Ea=(55.0506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_1H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['[CH2][C]=CCCF(6422)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(166.635,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC=[C]CF(11807)'],
    products = ['[CH2]CC=C([CH2])F(8231)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(98.9438,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #3972',
    isomers = [
        '[CH2]CC=[C]CF(11807)',
    ],
    reactants = [
        ('C2H4(30)', 'C#CCF(3587)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3972',
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

