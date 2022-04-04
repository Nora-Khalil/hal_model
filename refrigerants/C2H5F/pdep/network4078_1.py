species(
    label = 'FC[C]=CC=[C]CF(11137)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {6,S} {7,D} {13,S}
6  C u0 p0 c0 {5,S} {8,D} {14,S}
7  C u1 p0 c0 {3,S} {5,D}
8  C u1 p0 c0 {4,S} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (146.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([98,146,536,632,875,989,1009,1077,1202,1292,1337,1447,1384,1538,2913,3059,3013,3065,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.634514,'amu*angstrom^2'), symmetry=1, barrier=(14.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635398,'amu*angstrom^2'), symmetry=1, barrier=(14.6091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635655,'amu*angstrom^2'), symmetry=1, barrier=(14.615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.855961,0.0750828,-9.22772e-05,6.66962e-08,-2.029e-11,17734.1,26.4604], Tmin=(100,'K'), Tmax=(789.219,'K')), NASAPolynomial(coeffs=[8.63694,0.0356499,-1.73371e-05,3.39852e-09,-2.41009e-13,16505.8,-9.2372], Tmin=(789.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C#CCF(4379)',
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
    label = 'FCC1=CC=C1CF(11353)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u0 p0 c0 {5,D} {8,S} {13,S}
8  C u0 p0 c0 {6,D} {7,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-2.65783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817542,0.0655213,-5.43363e-05,2.26354e-08,-3.79974e-12,-201.598,21.9955], Tmin=(100,'K'), Tmax=(1406.34,'K')), NASAPolynomial(coeffs=[15.1152,0.0248542,-1.09597e-05,2.07246e-09,-1.44261e-13,-4222.96,-51.8572], Tmin=(1406.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.65783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'FC=C=CC=CCF(11140)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {8,D} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u0 p0 c0 {6,D} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-150.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864637,0.0660649,-5.73358e-05,2.56298e-08,-4.64741e-12,-18019.9,25.0732], Tmin=(100,'K'), Tmax=(1305.81,'K')), NASAPolynomial(coeffs=[13.9454,0.0259956,-1.13076e-05,2.13057e-09,-1.48426e-13,-21436.1,-41.524], Tmin=(1305.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCddFH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C=CC=[C]CF(11127)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {6,D} {10,S}
4  C u0 p0 c0 {3,S} {7,D} {11,S}
5  C u0 p0 c0 {7,D} {12,S} {13,S}
6  C u1 p0 c0 {2,S} {3,D}
7  C u0 p0 c0 {4,D} {5,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (260.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(0.980339,'amu*angstrom^2'), symmetry=1, barrier=(22.5399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983813,'amu*angstrom^2'), symmetry=1, barrier=(22.6198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22479,0.0605721,-5.5087e-05,2.65171e-08,-5.23266e-12,31469.2,22.7457], Tmin=(100,'K'), Tmax=(1199.72,'K')), NASAPolynomial(coeffs=[11.7583,0.0254523,-1.11769e-05,2.11686e-09,-1.48084e-13,28941.7,-29.9901], Tmin=(1199.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
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
    label = 'FC=C=CC=[C]CF(11363)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {7,D} {11,S}
5  C u0 p0 c0 {4,S} {8,D} {12,S}
6  C u0 p0 c0 {2,S} {8,D} {13,S}
7  C u1 p0 c0 {3,S} {4,D}
8  C u0 p0 c0 {5,D} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (87.0585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,2995,3025,975,1000,1300,1375,400,500,1630,1680,113,247,382,1207,3490,1685,370,540,610,2055,401.888,401.896],'cm^-1')),
        HinderedRotor(inertia=(0.120187,'amu*angstrom^2'), symmetry=1, barrier=(13.7748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120184,'amu*angstrom^2'), symmetry=1, barrier=(13.7749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05574,0.0679037,-7.20595e-05,4.15298e-08,-9.86526e-12,10574.2,24.8711], Tmin=(100,'K'), Tmax=(1005.55,'K')), NASAPolynomial(coeffs=[10.9538,0.0285297,-1.33239e-05,2.58853e-09,-1.83571e-13,8583.58,-22.9359], Tmin=(1005.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CF(4380)',
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
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,377.76,379.216,382.52,1900.36,1902.09,1902.61,1902.73,3346.48],'cm^-1')),
        HinderedRotor(inertia=(0.11674,'amu*angstrom^2'), symmetry=1, barrier=(12.0012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74483,0.0317272,-5.02124e-05,4.76765e-08,-1.76459e-11,37561.4,13.6835], Tmin=(100,'K'), Tmax=(834.526,'K')), NASAPolynomial(coeffs=[4.0464,0.0168537,-7.95764e-06,1.52224e-09,-1.05087e-13,37644.8,9.44116], Tmin=(834.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'FC[C]=CC#CCF(11364)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {6,D} {8,S} {13,S}
6  C u1 p0 c0 {3,S} {5,D}
7  C u0 p0 c0 {4,S} {8,T}
8  C u0 p0 c0 {5,S} {7,T}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (79.3416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,319,1023,1071,1259,1317,1409,3054,3019,3010,987.5,1337.5,450,1655,1685,370,2100,2250,500,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.902481,'amu*angstrom^2'), symmetry=1, barrier=(20.7498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.61913,'amu*angstrom^2'), symmetry=1, barrier=(60.219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900999,'amu*angstrom^2'), symmetry=1, barrier=(20.7157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34767,0.0617772,-5.88813e-05,3.15509e-08,-7.18609e-12,9635.24,25.4177], Tmin=(100,'K'), Tmax=(1028.14,'K')), NASAPolynomial(coeffs=[9.15467,0.0314036,-1.45675e-05,2.81675e-09,-1.99108e-13,8029.92,-12.4632], Tmin=(1028.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.3416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_S)"""),
)

species(
    label = 'FC=C[CH]C=[C]CF(11365)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {5,S} {6,S} {12,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {13,S}
7  C u0 p0 c0 {2,S} {5,D} {14,S}
8  C u1 p0 c0 {3,S} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (25.3007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15305,0.0521271,-1.75772e-05,-1.33699e-08,8.34268e-12,3154.53,28.947], Tmin=(100,'K'), Tmax=(1040.27,'K')), NASAPolynomial(coeffs=[13.7246,0.0254178,-1.02539e-05,1.92509e-09,-1.36681e-13,-631.414,-37.8253], Tmin=(1040.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.3007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + radical(C=CCJC=C) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
)

species(
    label = 'FC[C]=C[C]=CCF(11366)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,D} {13,S}
6  C u0 p0 c0 {7,D} {8,S} {14,S}
7  C u1 p0 c0 {4,S} {6,D}
8  C u1 p0 c0 {5,D} {6,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (107.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.911969,0.0752479,-8.66192e-05,4.92835e-08,-6.60094e-12,13058.5,26.1071], Tmin=(100,'K'), Tmax=(607.97,'K')), NASAPolynomial(coeffs=[8.65227,0.0356207,-1.67258e-05,3.20924e-09,-2.24258e-13,11908.5,-9.10088], Tmin=(607.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'FCC#C[CH][CH]CF(11367)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {6,S} {13,S}
6  C u1 p0 c0 {5,S} {8,S} {14,S}
7  C u0 p0 c0 {4,S} {8,T}
8  C u0 p0 c0 {6,S} {7,T}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (60.0367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.967557,0.0693239,-8.07051e-05,5.76153e-08,-1.7123e-11,7327.45,29.4549], Tmin=(100,'K'), Tmax=(843.394,'K')), NASAPolynomial(coeffs=[8.36191,0.0328437,-1.33154e-05,2.36349e-09,-1.57382e-13,6130.35,-4.66175], Tmin=(843.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.0367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(CsCsFHH) + group(CsCFHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cs_S) + radical(Sec_Propargyl)"""),
)

species(
    label = 'FC=[C]C=C[CH]CF(11368)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u1 p0 c0 {3,S} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {8,S} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (19.4706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.743681,0.0678803,-6.07037e-05,2.82973e-08,-5.33192e-12,2462.07,26.3145], Tmin=(100,'K'), Tmax=(1266.69,'K')), NASAPolynomial(coeffs=[14.121,0.025637,-1.06796e-05,1.96935e-09,-1.35714e-13,-926.909,-41.3858], Tmin=(1266.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.4706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Allyl_S) + radical(Cdj(Cd-CdH)(Cd-F1sH))"""),
)

species(
    label = '[CH2]C(F)=CC=[C]CF(11369)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {4,D} {7,S}
6  C u0 p0 c0 {4,S} {8,D} {12,S}
7  C u1 p0 c0 {5,S} {13,S} {14,S}
8  C u1 p0 c0 {3,S} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-6.23131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901636,0.0670206,-6.03313e-05,2.89714e-08,-5.71009e-12,-637.003,26.5072], Tmin=(100,'K'), Tmax=(1201.84,'K')), NASAPolynomial(coeffs=[12.4842,0.028471,-1.22179e-05,2.28246e-09,-1.58376e-13,-3421.06,-31.5013], Tmin=(1201.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.23131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C[C](F)CF(11370)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {8,S} {12,S}
7  C u0 p0 c0 {8,D} {13,S} {14,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-2.10095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.108,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.500727,0.0743622,-7.28271e-05,3.71583e-08,-7.61305e-12,-124.491,26.0432], Tmin=(100,'K'), Tmax=(1175.68,'K')), NASAPolynomial(coeffs=[14.7801,0.02578,-1.08435e-05,2.01095e-09,-1.39268e-13,-3482.12,-45.1579], Tmin=(1175.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.10095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CsCdCsF1s) + radical(C=CJC=C)"""),
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
    E0 = (14.8542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (23.1385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (23.2222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (264.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (171.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (226.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (167.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (492.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (154.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (211.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (165.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (164.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (102.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (146.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['C#CCF(4379)', 'C#CCF(4379)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['FCC1=CC=C1CF(11353)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleND_rad_out;CdsinglepriND_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['FC=C=CC=CCF(11140)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.284e+10,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'C=C=CC=[C]CF(11127)'],
    products = ['FC[C]=CC=[C]CF(11137)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.6964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(5)', 'FC=C=CC=[C]CF(11363)'],
    products = ['FC[C]=CC=[C]CF(11137)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.296016,'m^3/(mol*s)'), n=2.36777, Ea=(4.26157,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1084748591603159, var=0.5872085620069725, Tref=1000.0, N=10, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C#CCF(4379)', '[CH]=[C]CF(4380)'],
    products = ['FC[C]=CC=[C]CF(11137)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(539339,'m^3/(mol*s)'), n=0.286787, Ea=(38.5987,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.211724492539645, var=3.713057222477376, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Sp-2R!H#1R!H_N-Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Sp-2R!H#1R!H_N-Sp-4R!H-3R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'FC[C]=CC#CCF(11364)'],
    products = ['FC[C]=CC=[C]CF(11137)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1216.47,'m^3/(mol*s)'), n=1.464, Ea=(8.27115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_Ext-7R!H-R_N-Sp-7R!H-4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_Ext-7R!H-R_N-Sp-7R!H-4R!H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CF(4380)', '[CH]=[C]CF(4380)'],
    products = ['FC[C]=CC=[C]CF(11137)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.4099e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['FC=C[CH]C=[C]CF(11365)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.78195e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['FC[C]=C[C]=CCF(11366)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.5e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['FCC#C[CH][CH]CF(11367)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['FC=[C]C=C[CH]CF(11368)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.364e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R5HJ_3;Cd_rad_out_Cs;Cs_H_out_1H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['[CH2]C(F)=CC=[C]CF(11369)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.73244e+11,'s^-1'), n=0.48217, Ea=(88.0183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['FC[C]=CC=[C]CF(11137)'],
    products = ['C=[C]C=C[C](F)CF(11370)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(132.035,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #4078',
    isomers = [
        'FC[C]=CC=[C]CF(11137)',
    ],
    reactants = [
        ('C#CCF(4379)', 'C#CCF(4379)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4078',
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

