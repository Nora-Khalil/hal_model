species(
    label = 'F[C]=CC=[C]CF(11107)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {6,D} {10,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u1 p0 c0 {3,S} {4,D}
7  C u1 p0 c0 {2,S} {5,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (184.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.654183,'amu*angstrom^2'), symmetry=1, barrier=(15.041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652488,'amu*angstrom^2'), symmetry=1, barrier=(15.002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2307,0.0653649,-9.09006e-05,7.06601e-08,-2.23284e-11,22244.5,22.6531], Tmin=(100,'K'), Tmax=(771.391,'K')), NASAPolynomial(coeffs=[9.07373,0.0246865,-1.17825e-05,2.26804e-09,-1.58333e-13,21034.8,-13.1475], Tmin=(771.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'C2HF(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (95.331,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0062,'amu')),
        LinearRotor(inertia=(51.6236,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([429.793,429.793,596.357,596.357,1107.96,2365.05,3506.88],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1870.76,'J/mol'), sigma=(4.25,'angstroms'), dipoleMoment=(1,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4498,0.0030263,3.99146e-05,-8.9615e-08,5.74336e-11,11468.6,5.90915], Tmin=(10,'K'), Tmax=(555.749,'K')), NASAPolynomial(coeffs=[4.23833,0.0086714,-5.87678e-06,1.96876e-09,-2.53031e-13,11206.2,0.995297], Tmin=(555.749,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(95.331,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FCC1=CC=C1F(11249)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {6,S}
5  C u0 p0 c0 {4,D} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (1.79966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41438,0.0511863,-4.35308e-05,1.85016e-08,-3.12936e-12,314.189,20.3428], Tmin=(100,'K'), Tmax=(1414.5,'K')), NASAPolynomial(coeffs=[13.5808,0.0167807,-7.04462e-06,1.30495e-09,-8.99401e-14,-3127.61,-42.5719], Tmin=(1414.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.79966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F))"""),
)

species(
    label = 'FC=C=CC=CF(9015)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,D} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  C u0 p0 c0 {2,S} {7,D} {11,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-131.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14912,0.0548577,-5.06045e-05,2.34074e-08,-4.24413e-12,-15652.2,21.4569], Tmin=(100,'K'), Tmax=(1344.17,'K')), NASAPolynomial(coeffs=[14.802,0.014229,-5.26534e-06,9.20422e-10,-6.17951e-14,-19322.6,-48.4487], Tmin=(1344.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C=CC=[C]F(6371)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {5,D} {8,S}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u0 p0 c0 {5,D} {9,S} {10,S}
5  C u0 p0 c0 {2,D} {4,D}
6  C u1 p0 c0 {1,S} {3,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (298.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(1.14045,'amu*angstrom^2'), symmetry=1, barrier=(26.2212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60012,0.0507073,-5.25875e-05,2.82653e-08,-6.03222e-12,35979.8,18.253], Tmin=(100,'K'), Tmax=(1139.22,'K')), NASAPolynomial(coeffs=[11.6302,0.0154904,-6.2182e-06,1.13048e-09,-7.75765e-14,33694.5,-31.4437], Tmin=(1139.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cdj(Cd-CdH)(F1s))"""),
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
    label = 'F[C]=CC=C=CF(9016)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,S} {6,D} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {6,D} {10,S}
6  C u0 p0 c0 {3,D} {5,D}
7  C u1 p0 c0 {2,S} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (124.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,113,247,382,1207,3490,540,610,2055,167,640,1190,411.495,411.496],'cm^-1')),
        HinderedRotor(inertia=(0.111408,'amu*angstrom^2'), symmetry=1, barrier=(13.3868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42419,0.0581687,-7.01942e-05,4.427e-08,-1.11297e-11,15084.9,20.3992], Tmin=(100,'K'), Tmax=(969.471,'K')), NASAPolynomial(coeffs=[11.1621,0.0179907,-8.02997e-06,1.5225e-09,-1.06426e-13,13196.8,-26.2789], Tmin=(969.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = '[CH]=[C]F(1321)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u1 p0 c0 {3,D} {4,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (350.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([167,640,1190,1142.58,1502.03,3807.5],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.55591,0.0076799,-5.30098e-07,-4.54651e-09,2.16658e-12,42166.8,9.46363], Tmin=(100,'K'), Tmax=(1043.42,'K')), NASAPolynomial(coeffs=[5.86818,0.00351556,-1.29996e-06,2.62251e-10,-1.98915e-14,41428.4,-3.016], Tmin=(1043.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(CdCdF1s)"""),
)

species(
    label = 'F[C]=CC#CCF(11250)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {6,S} {7,D} {10,S}
5  C u0 p0 c0 {3,S} {6,T}
6  C u0 p0 c0 {4,S} {5,T}
7  C u1 p0 c0 {2,S} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (117.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([319,1023,1071,1259,1317,1409,3054,3019,3010,987.5,1337.5,450,1655,2100,2250,500,550,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.75478,'amu*angstrom^2'), symmetry=1, barrier=(40.3459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75174,'amu*angstrom^2'), symmetry=1, barrier=(40.2759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55959,0.0488336,-4.44661e-05,2.07017e-08,-3.834e-12,14273.2,22.3065], Tmin=(100,'K'), Tmax=(1303.27,'K')), NASAPolynomial(coeffs=[12.403,0.0155527,-6.16114e-06,1.10731e-09,-7.5278e-14,11446.8,-32.8789], Tmin=(1303.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCtH) + group(CdCFH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cdj(Cd-CtH)(F1s))"""),
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
    label = 'FC#CC=[C]CF(11251)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {5,D} {6,S} {10,S}
5  C u1 p0 c0 {3,S} {4,D}
6  C u0 p0 c0 {4,S} {7,T}
7  C u0 p0 c0 {2,S} {6,T}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (197.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([122,584,932,1043,1247,1392,1461,2986,3039,3010,987.5,1337.5,450,1655,1685,370,2175,525,239,401,1367,180],'cm^-1')),
        HinderedRotor(inertia=(0.983638,'amu*angstrom^2'), symmetry=1, barrier=(22.6158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0693096,'amu*angstrom^2'), symmetry=1, barrier=(12.47,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82981,0.0528872,-6.57955e-05,4.88605e-08,-1.55711e-11,23798.1,20.9448], Tmin=(100,'K'), Tmax=(748.505,'K')), NASAPolynomial(coeffs=[6.65287,0.027113,-1.41445e-05,2.85706e-09,-2.06068e-13,23076.1,-0.926581], Tmin=(748.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(CtCF) + radical(Cds_S)"""),
)

species(
    label = 'F[C]=C[CH]C=CF(9019)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  C u1 p0 c0 {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  C u1 p0 c0 {2,S} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (57.0109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,194,682,905,1196,1383,3221,167,640,1190,293.745,293.845],'cm^-1')),
        HinderedRotor(inertia=(0.00195619,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00195246,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69356,0.0413515,-1.0736e-05,-1.88462e-08,1.108e-11,6948.29,22.9717], Tmin=(100,'K'), Tmax=(969.539,'K')), NASAPolynomial(coeffs=[13.1648,0.0160199,-5.57349e-06,1.00262e-09,-7.1544e-14,3690.18,-37.3465], Tmin=(969.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.0109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(C=CCJC=C) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]C=C=CCF(11252)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {6,D} {7,S} {11,S}
6  C u0 p0 c0 {4,D} {5,D}
7  C u2 p0 c0 {2,S} {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (130.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07898,0.0511435,-3.35327e-05,-3.88312e-08,5.82024e-11,15725.5,19.7814], Tmin=(100,'K'), Tmax=(466.721,'K')), NASAPolynomial(coeffs=[6.12418,0.031082,-1.60041e-05,3.19065e-09,-2.27397e-13,15188.8,1.64379], Tmin=(466.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CsCCl_triplet)"""),
)

species(
    label = 'FC=[C]C=[C]CF(11253)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {6,D} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {7,D} {11,S}
6  C u1 p0 c0 {3,S} {4,D}
7  C u1 p0 c0 {4,S} {5,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (152.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34845,0.062338,-8.30195e-05,6.29162e-08,-1.95119e-11,18394.8,23.5278], Tmin=(100,'K'), Tmax=(783.859,'K')), NASAPolynomial(coeffs=[8.65647,0.0250455,-1.16562e-05,2.2223e-09,-1.54493e-13,17249.1,-9.94945], Tmin=(783.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(Cds_S) + radical(Cdj(Cd-CdH)(Cd-F1sH))"""),
)

species(
    label = 'F[C]=C=C[CH]CF(11254)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u1 p0 c0 {3,S} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {7,D}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (110.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78152,0.0499798,-4.3285e-05,1.99212e-08,-3.83263e-12,13348.1,20.8951], Tmin=(100,'K'), Tmax=(1207.96,'K')), NASAPolynomial(coeffs=[9.62285,0.0240143,-1.10421e-05,2.12655e-09,-1.49841e-13,11453.7,-18.4163], Tmin=(1207.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCsFHH) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(Cdj(Cdd-Cd)(F1s))"""),
)

species(
    label = 'F[CH][CH]C#CCF(11255)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
4  C u1 p0 c0 {5,S} {7,S} {10,S}
5  C u1 p0 c0 {2,S} {4,S} {11,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {4,S} {6,T}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (99.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69279,0.054182,-5.88368e-05,3.24291e-08,-4.91991e-12,12101.7,25.7309], Tmin=(100,'K'), Tmax=(674.272,'K')), NASAPolynomial(coeffs=[8.11087,0.024155,-9.93955e-06,1.78304e-09,-1.19802e-13,11053.3,-4.05956], Tmin=(674.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(CsCsFHH) + group(CsCFHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Csj(Cs-F1sHH)(Ct-Ct)(H)) + radical(Csj(Cs-CtHH)(F1s)(H))"""),
)

species(
    label = 'FC=[C][CH]C=CF(9023)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,D} {8,S}
4  C u1 p0 c0 {3,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  C u0 p0 c0 {2,S} {7,D} {11,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (44.7348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,615,860,1140,1343,3152,1685,370,484.737,485.034,1927.49],'cm^-1')),
        HinderedRotor(inertia=(0.674045,'amu*angstrom^2'), symmetry=1, barrier=(112.592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.673103,'amu*angstrom^2'), symmetry=1, barrier=(112.56,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67323,0.0431627,-1.87133e-05,-8.59964e-09,6.95782e-12,5471.23,23.0276], Tmin=(100,'K'), Tmax=(988.531,'K')), NASAPolynomial(coeffs=[12.4217,0.0173424,-6.34932e-06,1.14673e-09,-8.06558e-14,2482.74,-33.0708], Tmin=(988.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.7348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CdCFH) + group(CdCFH) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(F)=CC=[C]F(9069)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  C u0 p0 c0 {4,D} {5,S} {8,S}
4  C u0 p0 c0 {1,S} {3,D} {6,S}
5  C u0 p0 c0 {3,S} {7,D} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {5,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (31.3725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2783,0.0571416,-5.77883e-05,3.06711e-08,-6.49198e-12,3873.56,22.0096], Tmin=(100,'K'), Tmax=(1146.71,'K')), NASAPolynomial(coeffs=[12.3512,0.0185165,-7.26303e-06,1.29693e-09,-8.79342e-14,1334.1,-32.9264], Tmin=(1146.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.3725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(C=CC=CCJ) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = '[CH2][C]=CC=C(F)F(9070)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,D} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u1 p0 c0 {7,S} {10,S} {11,S}
7  C u1 p0 c0 {4,D} {6,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (20.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,182,240,577,636,1210,1413,3000,3100,440,815,1455,1000,1685,370,363.371],'cm^-1')),
        HinderedRotor(inertia=(1.03633,'amu*angstrom^2'), symmetry=1, barrier=(94.6856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01865,'amu*angstrom^2'), symmetry=1, barrier=(94.5136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44762,0.0496933,-3.66744e-05,9.79569e-09,4.84305e-13,2619.73,22.8554], Tmin=(100,'K'), Tmax=(1025.33,'K')), NASAPolynomial(coeffs=[12.6115,0.0176385,-6.60041e-06,1.17827e-09,-8.1154e-14,-273.985,-34.2306], Tmin=(1025.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFF) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    E0 = (19.1466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (27.4309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (44.1198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (268.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (175.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (199.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (172.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (269.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (250.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (497.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (158.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (215.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (181.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (169.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (204.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (52.1868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (107.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (144.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['C2HF(57)', 'C#CCF(4379)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['FCC1=CC=C1F(11249)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleND_rad_out;Cdsinglepri_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['FC=C=CC=CF(9015)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'C=C=CC=[C]F(6371)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.6964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(5)', 'F[C]=CC=C=CF(9016)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.296016,'m^3/(mol*s)'), n=2.36777, Ea=(4.26157,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1084748591603159, var=0.5872085620069725, Tref=1000.0, N=10, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]F(1321)', 'C#CCF(4379)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(6.03467,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'F[C]=CC#CCF(11250)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1216.47,'m^3/(mol*s)'), n=1.464, Ea=(8.23945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_Ext-7R!H-R_N-Sp-7R!H-4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H_Ext-4R!H-R_Ext-7R!H-R_N-Sp-7R!H-4R!H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2HF(57)', '[CH]=[C]CF(4380)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(539339,'m^3/(mol*s)'), n=0.286787, Ea=(26.9309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.211724492539645, var=3.713057222477376, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Sp-2R!H#1R!H_N-Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Sp-2R!H#1R!H_N-Sp-4R!H-3R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'FC#CC=[C]CF(11251)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(763,'m^3/(mol*s)'), n=1.64, Ea=(6.46828,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]F(1321)', '[CH]=[C]CF(4380)'],
    products = ['F[C]=CC=[C]CF(11107)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -37.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['F[C]=C[CH]C=CF(9019)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['F[C]C=C=CCF(11252)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['FC=[C]C=[C]CF(11253)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.50921e+08,'s^-1'), n=1.36125, Ea=(162.666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['F[C]=C=C[CH]CF(11254)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['F[CH][CH]C#CCF(11255)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(323312,'s^-1'), n=2.11016, Ea=(185.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_single;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['FC=[C][CH]C=CF(9023)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out] for rate rule [R5HJ_3;Cd_rad_out_single;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['[CH2]C(F)=CC=[C]F(9069)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(88.0183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C]=CC=[C]CF(11107)'],
    products = ['[CH2][C]=CC=C(F)F(9070)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(125.599,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #4014',
    isomers = [
        'F[C]=CC=[C]CF(11107)',
    ],
    reactants = [
        ('C2HF(57)', 'C#CCF(4379)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4014',
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

