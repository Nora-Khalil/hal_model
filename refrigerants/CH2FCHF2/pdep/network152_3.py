species(
    label = 'FC1=CC(F)=C1(307)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u0 p0 c0 {2,S} {3,D} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (29.0416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,180,1170.89,2106.02,2106.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89011,0.00686793,0.000107406,-2.40335e-07,1.58082e-10,3500.24,9.34061], Tmin=(10,'K'), Tmax=(520.912,'K')), NASAPolynomial(coeffs=[4.12106,0.0279909,-1.93509e-05,6.26873e-09,-7.66374e-13,3165.53,5.39524], Tmin=(520.912,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(29.0416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCC(F)DC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CC(F)[C]1(577)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {8,S}
5 C u0 p0 c0 {2,S} {4,D} {6,S}
6 C u0 p1 c0 {3,S} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (169.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,405,756,886,1212,205.073,205.191,205.197,205.222,205.26],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40673,0.0339548,-2.22855e-05,6.48333e-09,-7.36772e-13,20412.7,22.6086], Tmin=(100,'K'), Tmax=(2009.25,'K')), NASAPolynomial(coeffs=[13.0445,0.0127774,-6.47572e-06,1.23771e-09,-8.40933e-14,16137.9,-36.1352], Tmin=(2009.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(CdCCF) + group(CsJ2_singlet-CsH) + ring(Cyclobutene)"""),
)

species(
    label = '[C]=C(F)C=CF(863)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,D} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {8,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p1 c0 {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (65.3359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(1.49725,'amu*angstrom^2'), symmetry=1, barrier=(34.4248,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79287,0.051857,-7.20977e-05,5.19195e-08,-1.4898e-11,7934.65,15.5473], Tmin=(100,'K'), Tmax=(851.894,'K')), NASAPolynomial(coeffs=[9.58215,0.0152829,-7.69836e-06,1.52229e-09,-1.08225e-13,6607.52,-20.7829], Tmin=(851.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.3359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCCF) + group(CdJ2_singlet-Cds)"""),
)

species(
    label = 'FC1(F)[C]C=C1(862)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {4,D} {6,S} {8,S}
6 C u0 p1 c0 {3,S} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (144.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57806,0.033559,-2.33131e-05,7.67725e-09,-1.04052e-12,17417.8,22.7962], Tmin=(100,'K'), Tmax=(1616.93,'K')), NASAPolynomial(coeffs=[8.86336,0.0180102,-8.8888e-06,1.73008e-09,-1.21003e-13,15385.2,-10.547], Tmin=(1616.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + ring(Cyclobutene)"""),
)

species(
    label = 'FC1=CC=C1F(304)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {3,S} {6,D} {8,S}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (40.9583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89822,0.00612407,0.000102608,-2.16319e-07,1.34206e-10,4932.65,9.22518], Tmin=(10,'K'), Tmax=(550.566,'K')), NASAPolynomial(coeffs=[3.92734,0.0284965,-1.98739e-05,6.49561e-09,-8.00611e-13,4587.15,5.99356], Tmin=(550.566,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(40.9583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCCDC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=C=CC1F(864)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p0 c0 {3,S} {6,D} {8,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (128.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.64242,0.0242443,-1.21915e-05,2.13427e-09,-1.06039e-13,15487.1,11.0519], Tmin=(100,'K'), Tmax=(2632.43,'K')), NASAPolynomial(coeffs=[24.6872,-0.00137422,-1.21725e-06,2.727e-10,-1.63984e-14,2203.98,-115.031], Tmin=(2632.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(CdCddCF) + group(Cdd-CdsCds) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(cyclobutadiene_13)"""),
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
    label = 'FC1=C=C[C]1(865)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {4,S} {5,D} {6,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p1 c0 {2,S} {3,S}
5 C u0 p0 c0 {2,D} {3,D}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (741.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3010,987.5,1337.5,450,1655,324.887,325.216,325.315,325.346,325.389],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76327,0.0364618,-7.41586e-05,8.28743e-08,-3.35252e-11,89196.5,22.5416], Tmin=(100,'K'), Tmax=(845.604,'K')), NASAPolynomial(coeffs=[0.0963114,0.0257665,-1.38358e-05,2.71582e-09,-1.88402e-13,90481,39.8889], Tmin=(845.604,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCsH) + group(CdCddCF) + group(CsJ2_singlet-CsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'FC1[C]C#C1(866)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3 C u0 p1 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {5,T}
5 C u0 p0 c0 {3,S} {4,T}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (404.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,180,180,180,180,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43521,0.0370179,-5.25211e-05,4.25819e-08,-1.35904e-11,48696.6,4.47132], Tmin=(100,'K'), Tmax=(903.542,'K')), NASAPolynomial(coeffs=[6.00636,0.0151882,-6.28638e-06,1.09397e-09,-7.08363e-14,48297,-11.0353], Tmin=(903.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsJ2_singlet-CsH) + group(Ct-CtCs) + group(Ct-CtCs) + ring(Ring)"""),
)

species(
    label = 'FC1(F)C=C=C1(572)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {6,D} {7,S}
5 C u0 p0 c0 {3,S} {6,D} {8,S}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (77.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,516.465,518.581,519.329,519.781,519.894,520.574,523.013],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66204,0.0256297,-1.40514e-05,2.93739e-09,-2.21019e-13,9274.54,9.30355], Tmin=(100,'K'), Tmax=(3132.22,'K')), NASAPolynomial(coeffs=[32.239,-0.0088055,1.45336e-06,-1.52815e-10,8.8766e-15,-9637.3,-162.803], Tmin=(3132.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH]=C(F)C=[C]F(246)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u1 p0 c0 {2,S} {3,D}
6 C u1 p0 c0 {4,D} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (224.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,250,446,589,854,899,167,640,1190,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.960835,'amu*angstrom^2'), symmetry=1, barrier=(22.0915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.42,'J/mol'), sigma=(4.87902,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.03 K, Pc=58.54 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91288,0.0467648,-5.40569e-05,3.09687e-08,-6.99622e-12,27092,17.6804], Tmin=(100,'K'), Tmax=(1076.78,'K')), NASAPolynomial(coeffs=[11.1984,0.012271,-6.00549e-06,1.21852e-09,-8.89742e-14,25092.3,-27.8036], Tmin=(1076.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cdj(Cd-CdH)(F1s)) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = 'FC1=[C]C(F)[CH]1(573)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {2,S} {4,S} {6,D}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (174.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2950,1000,271,519,563,612,1379,1219.26,1219.26,1219.26],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32662,0.0295585,-9.28678e-06,-1.00841e-08,5.87345e-12,21010.5,17.1333], Tmin=(100,'K'), Tmax=(1050.45,'K')), NASAPolynomial(coeffs=[11.2427,0.0109863,-4.72753e-06,9.5983e-10,-7.19363e-14,18288.7,-30.3595], Tmin=(1050.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
)

species(
    label = 'F[C]1[C]=C(F)C1(574)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (184.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,246,474,533,1155,180,342.637,1425.04,1425.21,1425.34,1425.35],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25451,0.0336862,-2.53955e-05,9.12398e-09,-1.28831e-12,22296.9,16.3706], Tmin=(100,'K'), Tmax=(1677.76,'K')), NASAPolynomial(coeffs=[11.9705,0.0105218,-4.68529e-06,8.94644e-10,-6.20672e-14,19036.7,-35.5311], Tmin=(1677.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(Cd-CdH)_ring) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC1=C[C]=C1(575)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u0 p0 c0 {2,S} {5,D} {6,S}
4 C u0 p0 c0 {2,D} {5,S} {7,S}
5 C u1 p0 c0 {3,D} {4,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (438.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2750,3150,900,1100,180,382.604,1277.51,1278.01,1279.1,2494.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.70381,0.024069,-1.22472e-05,-4.58696e-10,1.46957e-12,52838.7,13.2228], Tmin=(100,'K'), Tmax=(1156.71,'K')), NASAPolynomial(coeffs=[8.81266,0.0106025,-4.71549e-06,9.24359e-10,-6.64605e-14,50913.1,-19.3529], Tmin=(1156.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclobutadiene-C1)"""),
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
    label = 'FC1=[C]C(F)=C1(576)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,D} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (243.169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,86,203,488,582,605,741,250,446,589,854,899,2541.68,2541.77],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40728,0.0290555,-1.5994e-05,-2.86546e-09,3.69721e-12,29309,15.2901], Tmin=(100,'K'), Tmax=(1023.72,'K')), NASAPolynomial(coeffs=[10.9888,0.00792502,-3.20232e-06,6.36893e-10,-4.77036e-14,26902.2,-29.4857], Tmin=(1023.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cdj(Cd-CdF1s)(Cd-CdF1s)_ring)"""),
)

species(
    label = 'FC1=C[C]C1F(578)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {2,S} {3,S} {5,D}
5 C u0 p0 c0 {4,D} {6,S} {8,S}
6 C u0 p1 c0 {3,S} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (178.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,2950,1000,342.411,345.117,348.58,348.809],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55226,0.0325161,-2.14989e-05,6.62543e-09,-8.246e-13,21519.4,23.5866], Tmin=(100,'K'), Tmax=(1780.58,'K')), NASAPolynomial(coeffs=[9.8236,0.0161812,-7.73787e-06,1.47313e-09,-1.01192e-13,18930,-15.6885], Tmin=(1780.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(Cyclobutene)"""),
)

species(
    label = 'FC1=C=C=C1(579)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u0 p0 c0 {1,S} {2,S} {5,D}
4 C u0 p0 c0 {2,D} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (469.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3010,987.5,1337.5,450,1655,485.723,485.748,485.802,485.851,485.865],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25774,0.0203048,-2.199e-06,-4.02444e-08,4.32603e-11,56488.8,9.56059], Tmin=(100,'K'), Tmax=(428.165,'K')), NASAPolynomial(coeffs=[4.34128,0.0161179,-8.32585e-06,1.67353e-09,-1.20341e-13,56341.6,4.61687], Tmin=(428.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCddCF) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[C]=C=CC(F)F(874)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {8,S}
5 C u0 p0 c0 {4,D} {6,D}
6 C u0 p1 c0 {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (110.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(0.456745,'amu*angstrom^2'), symmetry=1, barrier=(10.5015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11819,0.0464651,-7.40616e-05,6.83939e-08,-2.52559e-11,13347.7,16.8352], Tmin=(100,'K'), Tmax=(791.105,'K')), NASAPolynomial(coeffs=[5.38301,0.0216577,-1.12879e-05,2.23275e-09,-1.57254e-13,13090.9,3.49087], Tmin=(791.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(CdJ2_singlet-(Cdd-Cds))"""),
)

species(
    label = 'FC#CC=CF(875)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {4,D} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {8,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {2,S} {5,T}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-33.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83076,0.0135755,8.78156e-05,-2.57358e-07,2.14179e-10,-4026.92,10.7352], Tmin=(10,'K'), Tmax=(413.416,'K')), NASAPolynomial(coeffs=[4.06774,0.0270387,-1.8201e-05,5.79016e-09,-6.98588e-13,-4181.16,8.17277], Tmin=(413.416,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-33.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC#CCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[C]=C(F)C#C(876)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,D}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {6,S}
5 C u0 p1 c0 {2,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (433.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.65682,'amu*angstrom^2'), symmetry=1, barrier=(38.0935,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29971,0.0408226,-6.30138e-05,5.03325e-08,-1.60648e-11,52245.7,11.7643], Tmin=(100,'K'), Tmax=(766.748,'K')), NASAPolynomial(coeffs=[7.81799,0.0120342,-6.69386e-06,1.36308e-09,-9.79478e-14,51399.5,-13.3924], Tmin=(766.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCCF) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + group(CdJ2_singlet-Cds)"""),
)

species(
    label = '[C]=C=C=CF(877)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 C u0 p0 c0 {2,D} {4,D}
4 C u0 p0 c0 {3,D} {5,D}
5 C u0 p1 c0 {4,D}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (507.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,540,563.333,586.667,610,1970,2140,1302.11],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59322,0.0349483,-5.86619e-05,5.37798e-08,-1.91853e-11,61057.2,12.9311], Tmin=(100,'K'), Tmax=(831.361,'K')), NASAPolynomial(coeffs=[5.34236,0.014059,-7.14748e-06,1.38471e-09,-9.58123e-14,60864.9,1.7683], Tmin=(831.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCddFH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(CdJ2_singlet-(Cdd-Cds))"""),
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
    E0 = (324.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (46.4906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (188.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (113.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (227.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (461.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (368.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (220.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (91.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (60.7325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (106.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (370.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (313.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (188.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (238.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (256.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (96.1755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (333.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (276.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction20',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1(F)[C]C=C1(862)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(296.27,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]=C(F)C=CF(863)'],
    products = ['FC1=CC(F)[C]1(577)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1=CC=C1F(304)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(160.15,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction17',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.03468e+20,'s^-1'), n=-2.18552, Ea=(85.6267,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5822770821690705, var=6.058307656543806, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(199.312,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'FC1=C=C[C]1(865)'],
    products = ['FC1=CC(F)[C]1(577)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(142.189,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'FC1[C]C#C1(866)'],
    products = ['FC1=CC(F)[C]1(577)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(386.715,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC1(F)C=C=C1(572)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(284.614,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;Cdsingle_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['FC1=[C]C(F)[CH]1(573)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]1[C]=C(F)C1(574)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'FC1=C[C]=C1(575)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'FC1=[C]C(F)=C1(576)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.03232e+06,'m^3/(mol*s)'), n=-0.00929868, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Sp-6R!H-5R!H',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Sp-6R!H-5R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['FC1=C[C]C1F(578)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(151.739,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'FC1=C=C=C1(579)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(191.534,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]=C=CC(F)F(874)'],
    products = ['[C]=C(F)C=CF(863)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.45932e+11,'s^-1'), n=0.63878, Ea=(287.201,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]=C(F)C=CF(863)'],
    products = ['FC#CC=CF(875)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(172.098,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', '[C]=C(F)C#C(876)'],
    products = ['[C]=C(F)C=CF(863)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(322.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[C]=C=C=CF(877)'],
    products = ['[C]=C(F)C=CF(863)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(191.67,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #152',
    isomers = [
        'FC1=CC(F)=C1(307)',
        'FC1=CC(F)[C]1(577)',
        '[C]=C(F)C=CF(863)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #152',
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

