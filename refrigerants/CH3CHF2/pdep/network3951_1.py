species(
    label = '[O]CC=[C]CF(11800)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (105.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,122,584,932,1043,1247,1392,1461,2986,3039,3010,987.5,1337.5,450,1655,1685,370,180,1980.26],'cm^-1')),
        HinderedRotor(inertia=(0.160684,'amu*angstrom^2'), symmetry=1, barrier=(3.69443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160081,'amu*angstrom^2'), symmetry=1, barrier=(3.68058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27541,0.0436843,-5.65615e-05,5.71773e-08,-2.40845e-11,12722.3,23.3457], Tmin=(100,'K'), Tmax=(778.371,'K')), NASAPolynomial(coeffs=[1.33756,0.0360574,-1.78782e-05,3.50206e-09,-2.46712e-13,13245.4,30.0573], Tmin=(778.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = 'CH2O(20)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'FCC1=CCO1(11974)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-212.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08367,0.026683,3.31607e-05,-6.62056e-08,2.80398e-11,-25471.5,16.5767], Tmin=(100,'K'), Tmax=(963.549,'K')), NASAPolynomial(coeffs=[14.8302,0.0121367,-3.92477e-06,7.79954e-10,-6.25751e-14,-29709.1,-53.6872], Tmin=(963.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'O=CC=CCF(11980)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {6,S} {10,S}
6  C u0 p0 c0 {2,D} {5,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-292.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6911,0.0374559,-2.13957e-05,5.01853e-09,-4.40233e-13,-35185.6,16.5576], Tmin=(100,'K'), Tmax=(2668.66,'K')), NASAPolynomial(coeffs=[21.5489,0.00919058,-5.50858e-06,1.04976e-09,-6.84424e-14,-45250.7,-92.9311], Tmin=(2668.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-292.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'OCC=C=CF(7987)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {1,S} {6,D} {10,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-180.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10403,0.0440108,-3.88757e-05,2.00233e-08,-4.47826e-12,-21605.6,20.101], Tmin=(100,'K'), Tmax=(1034.28,'K')), NASAPolynomial(coeffs=[7.05911,0.0248468,-1.10816e-05,2.10741e-09,-1.47614e-13,-22630.5,-3.9713], Tmin=(1034.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
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
    label = 'O=CC=[C]CF(9609)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {6,D} {9,S}
5  C u0 p0 c0 {2,D} {4,S} {10,S}
6  C u1 p0 c0 {3,S} {4,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-55.0229,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.580684,'amu*angstrom^2'), symmetry=1, barrier=(13.3511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579897,'amu*angstrom^2'), symmetry=1, barrier=(13.333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0722,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0107,0.0495042,-7.25264e-05,7.09883e-08,-2.91072e-11,-6551.68,19.4992], Tmin=(100,'K'), Tmax=(731.819,'K')), NASAPolynomial(coeffs=[3.64949,0.0314754,-1.69791e-05,3.44784e-09,-2.4799e-13,-6548.62,13.7644], Tmin=(731.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.0229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
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
    label = 'C=C=CC[O](6612)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u0 p0 c0 {5,D} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {4,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (219.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.008603,'amu*angstrom^2'), symmetry=1, barrier=(9.86851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0818,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26302,0.0269559,-1.12238e-05,1.50775e-09,-7.90549e-15,26388.3,15.2474], Tmin=(100,'K'), Tmax=(2418.36,'K')), NASAPolynomial(coeffs=[17.8462,0.00937284,-4.37285e-06,7.36987e-10,-4.37843e-14,17423,-71.9387], Tmin=(2418.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = '[O]CC=C=CF(7988)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {1,S} {6,D} {10,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (45.5146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,113,247,382,1207,3490,540,610,2055,180,1480.55,1480.79],'cm^-1')),
        HinderedRotor(inertia=(0.0032195,'amu*angstrom^2'), symmetry=1, barrier=(5.01078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0722,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.68624,0.0286812,-1.3523e-05,2.45524e-09,-1.52838e-13,5461.89,15.133], Tmin=(100,'K'), Tmax=(3008.4,'K')), NASAPolynomial(coeffs=[35.7643,-0.00801279,1.80246e-06,-2.82667e-10,1.99845e-14,-16534.7,-179.437], Tmin=(3008.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.5146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = '[CH2][O](132)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88412,-0.00363928,3.2856e-05,-4.13634e-08,1.59641e-11,23210.8,7.4797], Tmin=(100,'K'), Tmax=(933.047,'K')), NASAPolynomial(coeffs=[6.69324,0.000290178,8.61305e-07,-1.56324e-10,7.33557e-15,21991.3,-9.60369], Tmin=(933.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[O]CC#CCF(11981)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,T}
6  C u0 p0 c0 {4,S} {5,T}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (26.9617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,319,1023,1071,1259,1317,1409,3054,3019,2100,2250,500,550,180,180,1340.6],'cm^-1')),
        HinderedRotor(inertia=(0.00544571,'amu*angstrom^2'), symmetry=1, barrier=(6.9318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.64358,'amu*angstrom^2'), symmetry=1, barrier=(60.7811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0722,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30863,0.0417786,-5.36196e-05,4.77332e-08,-1.71179e-11,3299.25,20.5388], Tmin=(100,'K'), Tmax=(875.062,'K')), NASAPolynomial(coeffs=[2.99445,0.028053,-1.19376e-05,2.14705e-09,-1.42882e-13,3584.7,19.6385], Tmin=(875.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.9617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtOsHH) + group(CsCFHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCOJ)"""),
)

species(
    label = 'O[CH]C=[C]CF(11982)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {6,D} {9,S}
5  C u1 p0 c0 {2,S} {4,S} {10,S}
6  C u1 p0 c0 {3,S} {4,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-3.10009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91619,0.0376525,-1.1323e-05,-1.13516e-08,6.51536e-12,-290.813,24.3782], Tmin=(100,'K'), Tmax=(1055.98,'K')), NASAPolynomial(coeffs=[11.5183,0.0181696,-7.6388e-06,1.46854e-09,-1.05584e-13,-3260.39,-26.9281], Tmin=(1055.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.10009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = '[O]C[CH]C=CF(7991)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {3,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {1,S} {5,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-25.5183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,220.803,221.51,2081.9],'cm^-1')),
        HinderedRotor(inertia=(1.44419,'amu*angstrom^2'), symmetry=1, barrier=(50.275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186916,'amu*angstrom^2'), symmetry=1, barrier=(50.2751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0089,0.0481185,-4.51582e-05,2.64124e-08,-7.0241e-12,-3001.16,18.8597], Tmin=(100,'K'), Tmax=(862.517,'K')), NASAPolynomial(coeffs=[5.78275,0.0306172,-1.47222e-05,2.88785e-09,-2.05629e-13,-3652.17,1.21109], Tmin=(862.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.5183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C[C]=CCF(11983)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (105.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,288.776,292.023,2291.52],'cm^-1')),
        HinderedRotor(inertia=(0.0869878,'amu*angstrom^2'), symmetry=1, barrier=(5.04829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858424,'amu*angstrom^2'), symmetry=1, barrier=(5.04647,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.11836,0.0289374,-1.16068e-05,1.13048e-09,6.338e-14,12592.9,15.1083], Tmin=(100,'K'), Tmax=(2698.4,'K')), NASAPolynomial(coeffs=[35.8083,-0.00513881,1.64853e-07,-6.20221e-12,4.55452e-15,-9205.92,-177.937], Tmin=(2698.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'OC[C]=[C]CF(11984)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {6,D}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (117.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,122,584,932,1043,1247,1392,1461,2986,3039,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.18247,'amu*angstrom^2'), symmetry=1, barrier=(4.19535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182572,'amu*angstrom^2'), symmetry=1, barrier=(4.19768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18246,'amu*angstrom^2'), symmetry=1, barrier=(4.19511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90171,0.0507231,-6.76597e-05,5.89202e-08,-2.15227e-11,14196.7,24.3033], Tmin=(100,'K'), Tmax=(778.062,'K')), NASAPolynomial(coeffs=[4.93617,0.0289734,-1.38736e-05,2.67646e-09,-1.86965e-13,13910.6,11.6215], Tmin=(778.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = '[O]C=C[CH]CF(11985)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {3,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {2,S} {5,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-131.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71555,0.037318,6.09235e-06,-4.04435e-08,1.95978e-11,-15730.8,19.3786], Tmin=(100,'K'), Tmax=(959.734,'K')), NASAPolynomial(coeffs=[15.3566,0.0120508,-3.7837e-06,7.08962e-10,-5.46523e-14,-19803.9,-53.4495], Tmin=(959.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = 'OC[CH][C]=CF(7994)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {3,S} {6,S} {9,S}
5  C u0 p0 c0 {1,S} {6,D} {10,S}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
"""),
    E0 = (-4.49022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57454,0.0519903,-4.80156e-05,2.31607e-08,-4.53537e-12,-451.623,21.505], Tmin=(100,'K'), Tmax=(1216.66,'K')), NASAPolynomial(coeffs=[11.2607,0.0201453,-8.75456e-06,1.64774e-09,-1.14882e-13,-2808.58,-27.1247], Tmin=(1216.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.49022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + radical(C=CCJCO) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = 'C=[C][CH]COF(8008)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {3,S} {6,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u1 p0 c0 {4,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (312.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,302.978,302.978,302.979],'cm^-1')),
        HinderedRotor(inertia=(0.376607,'amu*angstrom^2'), symmetry=1, barrier=(24.5324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376608,'amu*angstrom^2'), symmetry=1, barrier=(24.5324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894665,'amu*angstrom^2'), symmetry=1, barrier=(58.2787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5056,0.0571844,-6.00837e-05,3.39892e-08,-7.90687e-12,37723.1,19.587], Tmin=(100,'K'), Tmax=(1026.62,'K')), NASAPolynomial(coeffs=[10.1208,0.0236179,-1.10406e-05,2.14242e-09,-1.51815e-13,35954.2,-22.2031], Tmin=(1026.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=C(F)[CH]C[O](8009)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {3,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-42.2259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0802,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83229,0.0521405,-5.68388e-05,3.92394e-08,-1.19883e-11,-5004.56,18.7347], Tmin=(100,'K'), Tmax=(770.926,'K')), NASAPolynomial(coeffs=[5.99302,0.0305514,-1.48306e-05,2.9108e-09,-2.06916e-13,-5646.05,-0.255738], Tmin=(770.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.2259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(CdCsCdF) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CCJCO)"""),
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
    E0 = (7.73731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (327.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (16.0216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (71.1375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (32.7106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (105.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (77.0772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (257.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (163.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (190.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (166.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (407.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (123.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (147.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (241.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (113.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (158.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (65.4796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (324.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (98.6909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['CH2O(20)', 'C#CCF(3587)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[CH2]C=[C]CF(11864)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cd;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['FCC1=CCO1(11974)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['O=CC=CCF(11980)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['OCC=C=CF(7987)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2O(20)', '[CH]=[C]CF(4166)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.25713e-13,'m^3/(mol*s)'), n=5.31456, Ea=(9.82237,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.372338763425909, var=4.180287361148109, Tref=1000.0, N=384, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'O=CC=[C]CF(9609)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.5,'m^3/(mol*s)'), n=2.16, Ea=(17.8681,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_Sp-5R!H=4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0_4R!H->C_Sp-4C-1COS_Ext-4C-R_Sp-5R!H=4C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'C=C=CC[O](6612)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.5387,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]CC=C=CF(7988)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.296016,'m^3/(mol*s)'), n=2.36777, Ea=(4.20178,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1084748591603159, var=0.5872085620069725, Tref=1000.0, N=10, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][O](132)', 'C#CCF(3587)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.14341e+17,'m^3/(mol*s)'), n=-3.51414, Ea=(86.7854,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.599665927344668, var=38.10528332192646, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_N-1R!H->S_Sp-4R!H-3C_N-2R!H->S_Ext-2CNO-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]CC#CCF(11981)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6900,'m^3/(mol*s)'), n=1.64, Ea=(25.2655,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_N-Sp-7R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_N-Sp-7R!H=4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][O](132)', '[CH]=[C]CF(4166)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['O[CH]C=[C]CF(11982)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.33814e+11,'s^-1'), n=0.632482, Ea=(115.826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/OneDe] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['[O]C[CH]C=CF(7991)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C[C]=CCF(11983)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OC[C]=[C]CF(11984)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8074e+07,'s^-1'), n=1.23467, Ea=(93.7055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;O_H_out] + [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['[O]C=C[CH]CF(11985)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['OC[CH][C]=CF(7994)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.41845e+06,'s^-1'), n=1.36569, Ea=(57.7423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_1H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_1H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C][CH]COF(8008)'],
    products = ['[O]CC=[C]CF(11800)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(109.246,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]CC=[C]CF(11800)'],
    products = ['C=C(F)[CH]C[O](8009)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(90.9536,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #3951',
    isomers = [
        '[O]CC=[C]CF(11800)',
    ],
    reactants = [
        ('CH2O(20)', 'C#CCF(3587)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3951',
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

