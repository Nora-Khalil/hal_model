species(
    label = 'C=C([O])C1(F)[CH]O1(9320)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {2,S} {4,S} {8,S}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-129.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2950,1000,2950,3100,1380,975,1025,1650,394.855,394.944,395.106,395.121,395.356],'cm^-1')),
        HinderedRotor(inertia=(0.166849,'amu*angstrom^2'), symmetry=1, barrier=(18.5065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.983346,0.0578738,-6.17031e-05,3.16578e-08,-6.21393e-12,-15446,23.7194], Tmin=(100,'K'), Tmax=(1281.92,'K')), NASAPolynomial(coeffs=[16.6666,0.0081402,-2.57652e-06,4.23981e-10,-2.8152e-14,-19401.5,-55.5831], Tmin=(1281.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(CCsJO)"""),
)

species(
    label = 'FC1=CO1(3064)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-29.4903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,326,540,652,719,1357,788.854,788.882],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3010.76,'J/mol'), sigma=(4.86933,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=470.27 K, Pc=59.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.14414,0.0129866,6.60907e-06,-2.08037e-08,9.56931e-12,-3510.75,9.45209], Tmin=(100,'K'), Tmax=(957.696,'K')), NASAPolynomial(coeffs=[8.89006,0.00340757,-9.7292e-07,1.96224e-10,-1.667e-14,-5272.61,-21.4727], Tmin=(957.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.4903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + ring(Cyclopropene)"""),
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
    label = 'C=C([O])C1O[C]1F(9321)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-99.6197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,395,473,707,1436,2950,3100,1380,975,1025,1650,421.359,421.365,421.375,421.401,421.404,1606.36,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000949412,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3784.56,'J/mol'), sigma=(6.05541,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=591.14 K, Pc=38.67 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59507,0.045869,-3.41066e-05,5.9174e-09,2.46096e-12,-11888.6,25.4714], Tmin=(100,'K'), Tmax=(988.673,'K')), NASAPolynomial(coeffs=[13.8197,0.011203,-3.95551e-06,7.20187e-10,-5.1638e-14,-15028.9,-37.0225], Tmin=(988.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.6197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(CF)=C1[CH]O1(9543)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {7,S}
6  C u0 p0 c0 {3,S} {4,S} {5,D}
7  C u1 p0 c0 {2,S} {5,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-126.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,2950,1000,504.076,504.08,504.081,504.09,504.092,504.096,504.096,504.099,504.103],'cm^-1')),
        HinderedRotor(inertia=(0.000663399,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39991,0.0384637,1.71692e-05,-6.69477e-08,3.30916e-11,-15158.6,19.9046], Tmin=(100,'K'), Tmax=(926.827,'K')), NASAPolynomial(coeffs=[21.7617,-0.00112597,3.09193e-06,-6.08573e-10,3.46116e-14,-21007,-87.9706], Tmin=(926.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(C=C(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = '[O]C1=C(F)[CH]OC1(9544)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {1,S} {5,D} {7,S}
7  C u1 p0 c0 {2,S} {6,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-215.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14629,0.0359268,-1.95099e-05,4.10315e-09,-1.36321e-13,-25902.2,19.689], Tmin=(100,'K'), Tmax=(1509.29,'K')), NASAPolynomial(coeffs=[10.5603,0.0190791,-8.18385e-06,1.49351e-09,-1.00466e-13,-29063,-26.4245], Tmin=(1509.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(C=CCJ(O)C) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=[C]C1(F)[CH]O1(9491)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 C u0 p0 c0 {6,D} {8,S} {9,S}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (185.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,2950,3100,1380,975,1025,1650,1685,370,324.722,325.723,326.171,327.212],'cm^-1')),
        HinderedRotor(inertia=(0.312105,'amu*angstrom^2'), symmetry=1, barrier=(23.5393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94493,0.0435092,-4.18094e-05,2.02028e-08,-3.8969e-12,22350.6,19.1811], Tmin=(100,'K'), Tmax=(1244.98,'K')), NASAPolynomial(coeffs=[11.0975,0.0141027,-6.3791e-06,1.23028e-09,-8.70649e-14,20071.7,-26.98], Tmin=(1244.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(CCsJO) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OC2OC12F(9327)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {5,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-300.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70847,0.0329984,2.27379e-05,-5.97296e-08,2.59941e-11,-35993.7,13.191], Tmin=(100,'K'), Tmax=(991.621,'K')), NASAPolynomial(coeffs=[18.2046,0.00851419,-3.84506e-06,9.1366e-10,-7.79945e-14,-41333.1,-76.6808], Tmin=(991.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(CsCCFO) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC1([C]2CO2)[CH]O1(9545)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u1 p0 c0 {3,S} {4,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (71.1975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42228,0.0345514,-1.91904e-05,4.80415e-09,-4.6944e-13,8619.54,21.2805], Tmin=(100,'K'), Tmax=(2307.19,'K')), NASAPolynomial(coeffs=[13.6209,0.0151363,-6.56795e-06,1.15688e-09,-7.42361e-14,3452.04,-42.1089], Tmin=(2307.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.1975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Cs(F)(C)-Cs-O2s) + ring(Cs(C-F)-Cs-O2s) + radical(C2CsJO) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
)

species(
    label = '[O][C]1CC2OC12F(9546)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (21.4092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52909,0.0463677,-3.22904e-05,8.10014e-09,-1.02881e-13,2670.47,16.4904], Tmin=(100,'K'), Tmax=(1247.77,'K')), NASAPolynomial(coeffs=[14.3625,0.0155251,-7.59236e-06,1.51839e-09,-1.09349e-13,-1333.8,-51.4763], Tmin=(1247.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + polycyclic(s2_3_4_ane) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C([O])C(F)=C[O](9547)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u1 p2 c0 {5,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {3,S} {4,D} {8,S}
7  C u0 p0 c0 {5,D} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-217.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18206,'amu*angstrom^2'), symmetry=1, barrier=(27.1778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883647,0.0604524,-6.45058e-05,3.30566e-08,-6.49832e-12,-26075.5,22.3742], Tmin=(100,'K'), Tmax=(1257.6,'K')), NASAPolynomial(coeffs=[17.0429,0.00905551,-3.20269e-06,5.59416e-10,-3.82073e-14,-30140,-59.2889], Tmin=(1257.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C1([O])C2OC21F(9459)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (65.9309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17893,0.0602868,-6.25456e-05,3.11508e-08,-6.10482e-12,8032.67,16.5726], Tmin=(100,'K'), Tmax=(1233.59,'K')), NASAPolynomial(coeffs=[15.0451,0.0153248,-7.87348e-06,1.60451e-09,-1.16959e-13,4611.63,-53.2347], Tmin=(1233.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.9309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsHHH) + polycyclic(s2_3_3_ane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = 'C=C([O])C1=CO1(9548)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p0 c0 {2,S} {3,S} {6,D}
5 C u0 p0 c0 {1,S} {3,D} {7,S}
6 C u0 p0 c0 {4,D} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (126.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,1000,2950,3100,1380,975,1025,1650,384.341,384.357,384.361,384.364,384.39,384.391,384.392,384.399],'cm^-1')),
        HinderedRotor(inertia=(0.00114118,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0653,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.776756,0.0563966,-6.34337e-05,3.30818e-08,-6.32035e-12,15318.7,19.3042], Tmin=(100,'K'), Tmax=(1506.63,'K')), NASAPolynomial(coeffs=[17.2063,0.00195125,1.55045e-06,-4.4223e-10,3.37706e-14,11596.8,-62.6153], Tmin=(1506.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=O(338)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.3074e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'F[C]1[CH]O1(9215)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u1 p0 c0 {2,S} {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (97.7552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,661.567,661.569,2935.1],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18921,0.0134531,-5.29464e-07,-1.06374e-08,5.3494e-12,11790.3,13.5601], Tmin=(100,'K'), Tmax=(993.343,'K')), NASAPolynomial(coeffs=[8.19763,0.00359113,-1.19992e-06,2.57126e-10,-2.11232e-14,10286.8,-13.1285], Tmin=(993.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.7552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH]=C(O)C1(F)[CH]O1(9549)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {2,S} {4,S} {8,S}
7  C u1 p0 c0 {5,D} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-20.0973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2950,1000,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01749,'amu*angstrom^2'), symmetry=1, barrier=(23.3941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02779,'amu*angstrom^2'), symmetry=1, barrier=(23.6308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.368142,0.0674417,-7.81236e-05,4.14602e-08,-8.19286e-12,-2275.73,24.829], Tmin=(100,'K'), Tmax=(1375.29,'K')), NASAPolynomial(coeffs=[20.6007,0.00196178,5.29677e-07,-1.74118e-10,1.29992e-14,-7213.45,-76.9476], Tmin=(1375.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(CCsJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C1(F)CO1(9550)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u1 p0 c0 {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-66.0149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,350,440,435,1725,3120,650,792.5,1650,542.786,542.815,542.82,542.821,542.821],'cm^-1')),
        HinderedRotor(inertia=(0.0330212,'amu*angstrom^2'), symmetry=1, barrier=(6.90399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36395,0.0548839,-5.89766e-05,3.17439e-08,-6.69097e-12,-7842.18,22.5515], Tmin=(100,'K'), Tmax=(1160.61,'K')), NASAPolynomial(coeffs=[13.3778,0.0134786,-5.46318e-06,1.00516e-09,-6.96917e-14,-10630.8,-37.1976], Tmin=(1160.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.0149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(OF)=C1[CH]O1(9551)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {5,S}
4  C u0 p0 c0 {2,S} {5,D} {6,S}
5  C u0 p0 c0 {3,S} {4,D} {7,S}
6  C u1 p0 c0 {2,S} {4,S} {8,S}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (207.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,2950,1000,3000,3100,440,815,1455,1000,511.394,511.398,511.404,511.405,511.406,511.408,511.413,511.421],'cm^-1')),
        HinderedRotor(inertia=(0.000644561,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000644581,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22627,0.0433857,4.63829e-06,-5.58228e-08,2.9786e-11,25099.1,21.0459], Tmin=(100,'K'), Tmax=(920.044,'K')), NASAPolynomial(coeffs=[22.2627,-0.00204665,3.67028e-06,-7.47874e-10,4.59811e-14,19280.2,-89.2761], Tmin=(920.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])[C]1OC1F(9461)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u1 p0 c0 {2,S} {4,S} {6,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-114.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,350,440,435,1725,2950,3100,1380,975,1025,1650,504.063,504.063,504.066,504.067,504.067,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0365532,'amu*angstrom^2'), symmetry=1, barrier=(6.59053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67533,0.0440322,-3.04752e-05,3.71253e-09,2.7976e-12,-13629.4,24.8017], Tmin=(100,'K'), Tmax=(1004.19,'K')), NASAPolynomial(coeffs=[13.3034,0.0121872,-4.52658e-06,8.38583e-10,-6.01912e-14,-16694.5,-34.9793], Tmin=(1004.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(C2CsJO)"""),
)

species(
    label = 'O=[C]CC1(F)[CH]O1(9324)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6  C u1 p0 c0 {2,S} {4,S} {10,S}
7  C u1 p0 c0 {3,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-53.9898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,2850,1437.5,1250,1305,750,350,2950,1000,1855,455,950,180,180,1808.66,1809.24],'cm^-1')),
        HinderedRotor(inertia=(0.0767049,'amu*angstrom^2'), symmetry=1, barrier=(1.7636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0760707,'amu*angstrom^2'), symmetry=1, barrier=(1.74901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3643.24,'J/mol'), sigma=(6.15137,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.07 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89151,0.0497906,-5.80104e-05,3.95037e-08,-1.13041e-11,-6420.57,22.8758], Tmin=(100,'K'), Tmax=(837.559,'K')), NASAPolynomial(coeffs=[7.3804,0.0235769,-1.10637e-05,2.13582e-09,-1.50318e-13,-7340.02,-2.6319], Tmin=(837.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.9898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1=C(F)[CH]OO1(9552)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,D} {7,S}
5  C u0 p0 c0 {1,S} {4,D} {6,S}
6  C u1 p0 c0 {3,S} {5,S} {8,S}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (69.5198,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,2950,1000,3000,3100,440,815,1455,1000,180,999.166,1000.56,1000.69,1000.74,1000.98,1001.09,1001.1,1001.19,1001.71],'cm^-1')),
        HinderedRotor(inertia=(0.176174,'amu*angstrom^2'), symmetry=1, barrier=(4.05058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00396,0.0349688,-5.30101e-06,-1.72682e-08,8.43907e-12,8440.83,22.2938], Tmin=(100,'K'), Tmax=(1055.03,'K')), NASAPolynomial(coeffs=[11.7367,0.0173325,-7.61526e-06,1.50101e-09,-1.09551e-13,5315.05,-30.2634], Tmin=(1055.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.5198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + ring(12dioxolene) + radical(C=CCJO) + radical(C=C(O)CJ) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=[C]C1(F)[CH]O1(9538)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-76.5505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,1855,455,950,424.514,424.902,425.753],'cm^-1')),
        HinderedRotor(inertia=(0.163766,'amu*angstrom^2'), symmetry=1, barrier=(20.9821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0371,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12254,0.0361258,-3.36143e-05,1.21914e-08,-7.98239e-13,-9134.78,19.7479], Tmin=(100,'K'), Tmax=(1028.31,'K')), NASAPolynomial(coeffs=[12.3372,0.00549292,-2.20561e-06,4.35435e-10,-3.26256e-14,-11716.7,-32.1565], Tmin=(1028.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.5505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cs-Cs(F)-O2s) + radical(CCsJO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2OC12F(9380)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-309.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17344,0.0400931,-2.59125e-05,7.47862e-09,-8.46583e-13,-37190.8,13.9119], Tmin=(100,'K'), Tmax=(1999.01,'K')), NASAPolynomial(coeffs=[14.0019,0.0164247,-8.15267e-06,1.55581e-09,-1.05875e-13,-41920,-51.347], Tmin=(1999.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + longDistanceInteraction_cyclic(Cs(F)-CO) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = '[CH2][C]1OC2OC12F(9553)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {5,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
6  C u1 p0 c0 {3,S} {4,S} {7,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (8.63247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73058,0.0561976,-6.72367e-05,4.55134e-08,-1.31724e-11,1114.5,15.7222], Tmin=(100,'K'), Tmax=(818.449,'K')), NASAPolynomial(coeffs=[7.53739,0.0278173,-1.52215e-05,3.14331e-09,-2.29811e-13,164.012,-11.1287], Tmin=(818.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.63247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + polycyclic(s2_3_4_ane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
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
    E0 = (-89.6353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (60.9888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (181.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (59.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (428.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-120.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (102.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (22.0825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-35.3283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (67.8036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (228.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (131.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (66.3433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (258.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (172.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-21.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (302.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (68.1811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (192.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (371.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (305.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-120.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (10.0036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['FC1=CO1(3064)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(39.08,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 39.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(CF)=C1[CH]O1(9543)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(308.176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['[O]C1=C(F)[CH]OC1(9544)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(188.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'C=[C]C1(F)[CH]O1(9491)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['C=C1OC2OC12F(9327)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['FC1([C]2CO2)[CH]O1(9545)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['[O][C]1CC2OC12F(9546)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.3385e+08,'s^-1'), n=1.06125, Ea=(150.798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 149.5 to 150.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(F)=C[O](9547)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['[CH2]C1([O])C2OC21F(9459)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(196.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C=C([O])C1=CO1(9548)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(28.4415,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC1=CO1(3064)', '[CH2][C]=O(338)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00512641,'m^3/(mol*s)'), n=2.42326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7837337452902932, var=8.455757054931007, Tref=1000.0, N=46, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C]1[CH]O1(9215)', 'CH2CO(28)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.66944e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(28.7331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]1[CH]O1(9215)', '[CH2][C]=O(338)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(O)C1(F)[CH]O1(9549)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])C1(F)CO1(9550)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(OF)=C1[CH]O1(9551)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(94.5398,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])[C]1OC1F(9461)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(181.577,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]CC1(F)[CH]O1(9324)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C1=C(F)[CH]OO1(9552)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(301.531,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(T)(18)', 'O=[C]C1(F)[CH]O1(9538)'],
    products = ['C=C([O])C1(F)[CH]O1(9320)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['O=C1CC2OC12F(9380)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C1(F)[CH]O1(9320)'],
    products = ['[CH2][C]1OC2OC12F(9553)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.14898e+08,'s^-1'), n=0.935645, Ea=(138.719,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #3111',
    isomers = [
        'C=C([O])C1(F)[CH]O1(9320)',
    ],
    reactants = [
        ('FC1=CO1(3064)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3111',
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

