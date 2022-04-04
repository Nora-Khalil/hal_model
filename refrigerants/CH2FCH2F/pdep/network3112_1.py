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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59507,0.045869,-3.41066e-05,5.9174e-09,2.46096e-12,-11888.6,25.4714], Tmin=(100,'K'), Tmax=(988.673,'K')), NASAPolynomial(coeffs=[13.8197,0.011203,-3.95551e-06,7.20187e-10,-5.1638e-14,-15028.9,-37.0225], Tmin=(988.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.6197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(CsCsF1sO2s)"""),
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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.983346,0.0578738,-6.17031e-05,3.16578e-08,-6.21393e-12,-15446,23.7194], Tmin=(100,'K'), Tmax=(1281.92,'K')), NASAPolynomial(coeffs=[16.6666,0.0081402,-2.57652e-06,4.23981e-10,-2.8152e-14,-19401.5,-55.5831], Tmin=(1281.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(CCsJO)"""),
)

species(
    label = '[O]C1=C[C](F)OC1(9453)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u1 p2 c0 {5,S}
4  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {7,S} {10,S}
7  C u1 p0 c0 {1,S} {2,S} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-183.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88483,0.0344805,5.53444e-06,-3.65684e-08,1.74832e-11,-21982.2,19.1473], Tmin=(100,'K'), Tmax=(970.932,'K')), NASAPolynomial(coeffs=[14.4698,0.011653,-4.03149e-06,7.82787e-10,-6.03304e-14,-25793.9,-48.2405], Tmin=(970.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(Cs_P)"""),
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
    label = 'C=[C]C1O[C]1F(9454)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 C u0 p0 c0 {6,D} {8,S} {9,S}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (214.976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,2950,3100,1380,975,1025,1650,1685,370,180,348.746,1081.67,1081.75,1081.84,4000],'cm^-1')),
        HinderedRotor(inertia=(0.577695,'amu*angstrom^2'), symmetry=1, barrier=(13.2823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0643,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20468,0.0355073,-2.74894e-05,1.06085e-08,-1.64201e-12,25923.5,22.2053], Tmin=(100,'K'), Tmax=(1522.77,'K')), NASAPolynomial(coeffs=[10.6506,0.0133216,-5.63533e-06,1.04084e-09,-7.12348e-14,23351.2,-22.0932], Tmin=(1522.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OC2(F)OC12(9326)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {5,S} {6,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-289.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73019,0.0291706,4.00137e-05,-8.29939e-08,3.57158e-11,-34702.6,14.0478], Tmin=(100,'K'), Tmax=(967.43,'K')), NASAPolynomial(coeffs=[20.202,0.0042663,-1.17696e-06,3.85464e-10,-4.22695e-14,-40685.3,-86.9047], Tmin=(967.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFOO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'C=C(O)C1=C(F)O1(9455)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {10,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {1,S} {2,S} {4,D}
7  C u0 p0 c0 {5,D} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-196.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785967,0.0673313,-8.35447e-05,4.99185e-08,-1.14369e-11,-23540.3,19.1788], Tmin=(100,'K'), Tmax=(1081.42,'K')), NASAPolynomial(coeffs=[16.368,0.00969519,-3.5986e-06,6.33185e-10,-4.30814e-14,-26910.4,-57.2148], Tmin=(1081.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFO) + group(Cds-CdsHH) + ring(Cyclopropene)"""),
)

species(
    label = 'F[C]1OC1[C]1CO1(9456)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u1 p0 c0 {1,S} {3,S} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (81.0143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54393,0.0463958,-4.3345e-05,2.27531e-08,-4.62399e-12,9838.73,22.5839], Tmin=(100,'K'), Tmax=(1398.2,'K')), NASAPolynomial(coeffs=[9.93991,0.0157774,-3.41791e-06,3.40295e-10,-1.30069e-14,8135.92,-18.4293], Tmin=(1398.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.0143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(CsCFHO) + ring(Ethylene_oxide) + ring(O2s-Cs(F)-Cs(C)) + radical(C2CsJO) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = '[O][C]1CC2(F)OC12(9457)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (21.4092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52909,0.0463677,-3.22904e-05,8.10014e-09,-1.02881e-13,2670.47,16.4904], Tmin=(100,'K'), Tmax=(1247.77,'K')), NASAPolynomial(coeffs=[14.3625,0.0155251,-7.59236e-06,1.51839e-09,-1.09349e-13,-1333.8,-51.4763], Tmin=(1247.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + polycyclic(s2_3_4_ane) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C([O])=CC(=O)F(9458)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,D} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {4,D} {6,S}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-273.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,255,533,799,832,1228,342.103,2507.88],'cm^-1')),
        HinderedRotor(inertia=(0.140879,'amu*angstrom^2'), symmetry=1, barrier=(11.7112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140835,'amu*angstrom^2'), symmetry=1, barrier=(11.7165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36951,0.057433,-6.93855e-05,4.31564e-08,-1.05421e-11,-32849.5,23.3105], Tmin=(100,'K'), Tmax=(1005.23,'K')), NASAPolynomial(coeffs=[12.061,0.0148893,-5.90179e-06,1.054e-09,-7.12675e-14,-34999,-28.3257], Tmin=(1005.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(COCFO) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
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
    label = 'C=C([O])C1=C(F)O1(9460)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {2,S} {5,S} {6,D}
5 C u0 p0 c0 {3,S} {4,S} {7,D}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 C u0 p0 c0 {5,D} {8,S} {9,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-58.9068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,326,540,652,719,1357,2950,3100,1380,975,1025,1650,378.982,379.222,379.571,379.753,380.075],'cm^-1')),
        HinderedRotor(inertia=(0.0683739,'amu*angstrom^2'), symmetry=1, barrier=(7.02963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27021,0.0619014,-8.57057e-05,5.97697e-08,-1.62823e-11,-6988.18,18.9468], Tmin=(100,'K'), Tmax=(904.439,'K')), NASAPolynomial(coeffs=[12.1124,0.0139503,-6.17919e-06,1.15015e-09,-7.89834e-14,-8949.39,-32.2715], Tmin=(904.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.9068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFO) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=C(C)OJ)"""),
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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67533,0.0440322,-3.04752e-05,3.71253e-09,2.7976e-12,-13629.4,24.8017], Tmin=(100,'K'), Tmax=(1004.19,'K')), NASAPolynomial(coeffs=[13.3034,0.0121872,-4.52658e-06,8.38583e-10,-6.01912e-14,-16694.5,-34.9793], Tmin=(1004.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(C2CsJO)"""),
)

species(
    label = 'C=C(O)[C]1O[C]1F(9462)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {10,S}
4  C u1 p0 c0 {2,S} {5,S} {6,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u0 p0 c0 {5,D} {8,S} {9,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-57.3197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,395,473,707,1436,2950,3100,1380,975,1025,1650,440.537,440.623,440.698,440.735,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0859574,'amu*angstrom^2'), symmetry=1, barrier=(11.8419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859363,'amu*angstrom^2'), symmetry=1, barrier=(11.8415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17687,0.0564881,-6.30603e-05,3.42117e-08,-7.13259e-12,-6787.53,26.5863], Tmin=(100,'K'), Tmax=(1186.24,'K')), NASAPolynomial(coeffs=[15.2274,0.00910947,-3.14977e-06,5.41879e-10,-3.66431e-14,-10121,-43.5992], Tmin=(1186.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.3197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C2CsJO) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH]=C(O)C1O[C]1F(9463)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {9,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u1 p0 c0 {5,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (9.67169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,1000,350,440,435,1725,395,473,707,1436,3120,650,792.5,1650,347.611,348.714,348.981,352.079,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00136868,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241027,'amu*angstrom^2'), symmetry=1, barrier=(20.815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635192,0.0593383,-6.34274e-05,3.14187e-08,-5.78076e-12,1296.94,27.8287], Tmin=(100,'K'), Tmax=(1514.83,'K')), NASAPolynomial(coeffs=[18.7937,0.0032497,1.72385e-07,-1.18409e-10,9.38487e-15,-3270.53,-64.2342], Tmin=(1514.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.67169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C1OC1F(9464)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u1 p0 c0 {6,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-47.0776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,261,493,600,1152,1365,1422,3097,350,440,435,1725,3120,650,792.5,1650,507.339,508.26,509.225,509.522,510.712,3759.61],'cm^-1')),
        HinderedRotor(inertia=(0.000651986,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61003,0.0414087,-1.23192e-05,-2.23122e-08,1.37519e-11,-5565.94,24.3249], Tmin=(100,'K'), Tmax=(953.785,'K')), NASAPolynomial(coeffs=[16.148,0.00763972,-1.98902e-06,3.67607e-10,-3.00763e-14,-9576.37,-51.6102], Tmin=(953.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.0776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'O=[C]CC1O[C]1F(9325)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {2,S} {4,S}
7  C u1 p0 c0 {3,D} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-39.0225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,1855,455,950,180,180,180,1226.53,1233.1,2741.74],'cm^-1')),
        HinderedRotor(inertia=(0.00293732,'amu*angstrom^2'), symmetry=1, barrier=(15.7346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289572,'amu*angstrom^2'), symmetry=1, barrier=(46.9664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3566.27,'J/mol'), sigma=(5.7756,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.04 K, Pc=42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10311,0.0444462,-4.19156e-05,2.23718e-08,-5.11009e-12,-4627.29,22.9301], Tmin=(100,'K'), Tmax=(1021.15,'K')), NASAPolynomial(coeffs=[7.4971,0.0233171,-1.08784e-05,2.10885e-09,-1.4927e-13,-5728.9,-3.20565], Tmin=(1021.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.0225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCFHO) + group(Cds-OdCsH) + ring(O2s-Cs(F)-Cs(C)) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1[CH][C](F)OO1(9465)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,S} {7,D}
5  C u1 p0 c0 {4,S} {6,S} {8,S}
6  C u1 p0 c0 {1,S} {3,S} {5,S}
7  C u0 p0 c0 {4,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (65.5235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,2950,3100,1380,975,1025,1650,637.861,637.863,637.863,637.865,637.865,637.866,637.866,637.866,637.867,637.867,637.867,637.868],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91152,0.0325331,1.50489e-05,-4.45051e-08,1.9162e-11,7967.68,19.8986], Tmin=(100,'K'), Tmax=(1004.17,'K')), NASAPolynomial(coeffs=[14.2628,0.015139,-6.47848e-06,1.3289e-09,-1.01668e-13,3883.56,-47.7249], Tmin=(1004.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJCO) + radical(CsCsF1sO2s)"""),
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
    label = 'O=[C]C1O[C]1F(9466)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-46.7816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,1855,455,950,605.14,605.616,605.684,3681.97,4000],'cm^-1')),
        HinderedRotor(inertia=(0.109413,'amu*angstrom^2'), symmetry=1, barrier=(28.5146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0371,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55131,0.0262177,-1.30506e-05,-4.88443e-09,4.38553e-12,-5569.37,22.1602], Tmin=(100,'K'), Tmax=(1000.83,'K')), NASAPolynomial(coeffs=[10.4433,0.00699753,-2.71134e-06,5.29695e-10,-3.96296e-14,-7766.19,-19.0038], Tmin=(1000.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.7816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsH) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2(F)OC12(9390)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-326.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14107,0.0288498,1.54273e-05,-3.87993e-08,1.54843e-11,-39137.1,15.6159], Tmin=(100,'K'), Tmax=(1057.8,'K')), NASAPolynomial(coeffs=[12.7611,0.0172189,-8.53372e-06,1.79785e-09,-1.36202e-13,-42979.9,-43.7603], Tmin=(1057.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'CC(=O)C1=C(F)O1(9467)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {6,S} {7,D}
6  C u0 p0 c0 {3,D} {4,S} {5,S}
7  C u0 p0 c0 {1,S} {2,S} {5,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-195.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90909,0.0510767,-6.36466e-05,5.00449e-08,-1.72846e-11,-23490.5,18.7553], Tmin=(100,'K'), Tmax=(688.202,'K')), NASAPolynomial(coeffs=[5.73842,0.0288222,-1.51468e-05,3.0682e-09,-2.21596e-13,-24017.6,1.7115], Tmin=(688.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCFO) + ring(Cyclopropene)"""),
)

species(
    label = '[CH2][C]1OC2(F)OC12(9468)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {5,S} {6,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6  C u1 p0 c0 {3,S} {4,S} {7,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (18.5484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4645,0.0614981,-8.50589e-05,6.68376e-08,-2.18465e-11,2316.89,16.0254], Tmin=(100,'K'), Tmax=(739.006,'K')), NASAPolynomial(coeffs=[7.92061,0.0265557,-1.41393e-05,2.86452e-09,-2.06353e-13,1362.61,-13.1694], Tmin=(739.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.5484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCFOO) + group(Cs-CsHHH) + polycyclic(s2_3_4_ane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = 'CC(=O)[C]1O[C]1F(9469)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,D} {4,S} {5,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-81.5082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,395,473,707,1436,195.935,1175.45,1177.94,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(1.00208,'amu*angstrom^2'), symmetry=1, barrier=(27.6585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199345,'amu*angstrom^2'), symmetry=1, barrier=(5.23589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29245,0.0407002,-3.50812e-05,1.66978e-08,-3.44484e-12,-9744.34,24.982], Tmin=(100,'K'), Tmax=(1104.39,'K')), NASAPolynomial(coeffs=[7.18797,0.022969,-1.09982e-05,2.15997e-09,-1.53911e-13,-10825.6,0.877932], Tmin=(1104.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.5082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + ring(Cs-Cs(F)-O2s) + radical(C2CsJO) + radical(CsCsF1sO2s)"""),
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
    E0 = (-78.7679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (71.8561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (101.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (469.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-79.7947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-24.6789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (143.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (37.0852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-75.8405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (80.1423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (142.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (164.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (84.5515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (269.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (69.6873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (78.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (213.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (8.77158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (218.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (392.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (346.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-79.7947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-24.6789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (154.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (113.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['FC1=CO1(3064)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(9.31107,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 9.3 kJ/mol to match endothermicity of reaction."""),
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
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['[O]C1=C[C](F)OC1(9453)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(190.068,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'C=[C]C1O[C]1F(9454)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['C=C1OC2(F)OC12(9326)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['C=C(O)C1=C(F)O1(9455)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['F[C]1OC1[C]1CO1(9456)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['[O][C]1CC2(F)OC12(9457)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([O])=CC(=O)F(9458)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['[CH2]C1([O])C2OC21F(9459)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(168.221,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['FC1=CO1(3064)', '[CH2][C]=O(338)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00512641,'m^3/(mol*s)'), n=2.42326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7837337452902932, var=8.455757054931007, Tref=1000.0, N=46, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_5R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'C=C([O])C1=C(F)O1(9460)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C]1[CH]O1(9215)', 'CH2CO(28)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.66944e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(36.074,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]1[CH]O1(9215)', '[CH2][C]=O(338)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['C=C([O])[C]1OC1F(9461)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.71585e+08,'s^-1'), n=1.23602, Ea=(157.766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_noH;Cs_H_out_noH] + [R2H_S_cy3;C_rad_out_single;Cs_H_out_noH] + [R2H_S_cy3;C_rad_out_noH;Cs_H_out] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)[C]1O[C]1F(9462)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.69782e+07,'s^-1'), n=1.4084, Ea=(124.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C1O[C]1F(9463)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C([O])C1OC1F(9464)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]CC1O[C]1F(9325)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C1[CH][C](F)OO1(9465)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(315.785,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(T)(18)', 'O=[C]C1O[C]1F(9466)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['O=C1CC2(F)OC12(9390)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['CC(=O)C1=C(F)O1(9467)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C1O[C]1F(9321)'],
    products = ['[CH2][C]1OC2(F)OC12(9468)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(=O)[C]1O[C]1F(9469)'],
    products = ['C=C([O])C1O[C]1F(9321)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.22635e+07,'s^-1'), n=1.82798, Ea=(183.706,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3112',
    isomers = [
        'C=C([O])C1O[C]1F(9321)',
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
    label = 'PDepNetwork #3112',
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

