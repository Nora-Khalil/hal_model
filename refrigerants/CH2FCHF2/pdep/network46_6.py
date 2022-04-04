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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79287,0.051857,-7.20977e-05,5.19195e-08,-1.4898e-11,7934.65,15.5473], Tmin=(100,'K'), Tmax=(851.894,'K')), NASAPolynomial(coeffs=[9.58215,0.0152829,-7.69836e-06,1.52229e-09,-1.08225e-13,6607.52,-20.7829], Tmin=(851.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.3359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(CdCCF) + group(CdJ2_singlet-Cds)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55226,0.0325161,-2.14989e-05,6.62543e-09,-8.246e-13,21519.4,23.5866], Tmin=(100,'K'), Tmax=(1780.58,'K')), NASAPolynomial(coeffs=[9.8236,0.0161812,-7.73787e-06,1.47313e-09,-1.01192e-13,18930,-15.6885], Tmin=(1780.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(Cyclobutene)"""),
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
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,145,326,398,834,1303,180,771.73,1631.86,1636.88],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3176.46,'J/mol'), sigma=(4.96469,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=496.16 K, Pc=58.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.64242,0.0242443,-1.21915e-05,2.13427e-09,-1.06039e-13,15487.1,11.0519], Tmin=(100,'K'), Tmax=(2632.43,'K')), NASAPolynomial(coeffs=[24.6872,-0.00137422,-1.21725e-06,2.727e-10,-1.63984e-14,2203.98,-115.031], Tmin=(2632.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(CdCddCF) + group(Cdd-CdsCds) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'C2HF(58)',
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
    label = '[CH]=[C]F(252)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.55592,0.00767983,-5.29841e-07,-4.54684e-09,2.16672e-12,42166.8,9.46361], Tmin=(100,'K'), Tmax=(1043.41,'K')), NASAPolynomial(coeffs=[5.86815,0.00351561,-1.29999e-06,2.62257e-10,-1.9892e-14,41428.4,-3.01582], Tmin=(1043.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(CdCdF1s)"""),
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
    label = '[CH]=C(F)C#CF(308)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p0 c0 {3,S} {6,T}
5 C u1 p0 c0 {3,D} {7,S}
6 C u0 p0 c0 {2,S} {4,T}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (227.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([322,421,608,607,2175,525,3120,650,792.5,1650,239,401,1367,334.182],'cm^-1')),
        HinderedRotor(inertia=(0.0015111,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83565,0.0113127,0.000111388,-3.17083e-07,2.50404e-10,27353.3,11.2409], Tmin=(10,'K'), Tmax=(459.087,'K')), NASAPolynomial(coeffs=[6.13844,0.0205514,-1.45402e-05,4.81769e-09,-6.01006e-13,26833,-1.43907], Tmin=(459.087,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(227.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""[CH]DC(F)C#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C#CC=[C]F(309)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u0 p0 c0 {2,S} {5,T}
4 C u1 p0 c0 {1,S} {2,D}
5 C u0 p0 c0 {3,T} {7,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (341.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2175,525,167,640,1190,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.55519,'amu*angstrom^2'), symmetry=1, barrier=(35.7568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89641,0.00672961,9.14697e-05,-2.2232e-07,1.57074e-10,41055.9,10.0413], Tmin=(10,'K'), Tmax=(493.725,'K')), NASAPolynomial(coeffs=[4.70259,0.0208007,-1.3873e-05,4.43968e-09,-5.42535e-13,40725.2,4.17794], Tmin=(493.725,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(341.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""C#CCD[C]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]=C=C=[C]F(310)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 C u0 p0 c0 {3,D} {4,D}
3 C u0 p0 c0 {2,D} {5,D}
4 C u1 p0 c0 {2,D} {6,S}
5 C u1 p0 c0 {1,S} {3,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (509.626,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,137,207,812,180,1590.06,2099.35],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40428,0.0378747,-6.11711e-05,5.05318e-08,-1.59243e-11,61348.6,15.2978], Tmin=(100,'K'), Tmax=(903.844,'K')), NASAPolynomial(coeffs=[7.43498,0.00944976,-3.77241e-06,6.52976e-10,-4.18301e-14,60690.9,-7.07156], Tmin=(903.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=C=CJ) + radical(Cdj(Cdd-Cdd)(F1s))"""),
)

species(
    label = '[CH]C(F)=C=CF(311)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {5,D} {6,S}
4 C u0 p0 c0 {2,S} {5,D} {7,S}
5 C u0 p0 c0 {3,D} {4,D}
6 C u2 p0 c0 {3,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (149.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32587,0.0400923,-3.38632e-05,1.64559e-08,-3.54973e-12,18091.8,17.0495], Tmin=(100,'K'), Tmax=(1047.8,'K')), NASAPolynomial(coeffs=[6.41226,0.0244923,-1.15305e-05,2.24654e-09,-1.59417e-13,17235.5,-2.85571], Tmin=(1047.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(F)=C=[C]F(293)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 C u0 p0 c0 {3,D} {6,D}
6 C u1 p0 c0 {2,S} {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (138.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([132,421,945,1169,1268,3000,3100,440,815,1455,1000,540,610,2055,137,207,812],'cm^-1')),
        HinderedRotor(inertia=(0.176613,'amu*angstrom^2'), symmetry=1, barrier=(109.127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25385,0.0399254,-4.122e-05,2.24796e-08,-4.99692e-12,16720.6,17.1289], Tmin=(100,'K'), Tmax=(1074.93,'K')), NASAPolynomial(coeffs=[8.78759,0.0156121,-7.29215e-06,1.43763e-09,-1.0309e-13,15315.9,-14.8646], Tmin=(1074.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Allyl_P) + radical(Cdj(Cdd-Cd)(F1s))"""),
)

species(
    label = '[CH]=C=C[C](F)F(282)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 C u0 p0 c0 {3,D} {6,D}
6 C u1 p0 c0 {5,D} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (63.4463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([161,297,490,584,780,1358,540,610,2055,3120,650,792.5,1650,180,1376.15,1376.75,1376.83],'cm^-1')),
        HinderedRotor(inertia=(0.275045,'amu*angstrom^2'), symmetry=1, barrier=(6.32383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94371,0.0481008,-6.95106e-05,5.52628e-08,-1.75032e-11,7702.14,18.9919], Tmin=(100,'K'), Tmax=(822.058,'K')), NASAPolynomial(coeffs=[8.05892,0.0161101,-7.05893e-06,1.3086e-09,-8.90779e-14,6772.25,-8.85281], Tmin=(822.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.4463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCdF1sF1s) + radical(C=C=CJ)"""),
)

species(
    label = 'F[C]C=C=CF(306)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {5,D} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {5,D} {8,S}
5 C u0 p0 c0 {3,D} {4,D}
6 C u2 p0 c0 {2,S} {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (168.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3189,0.0411507,-5.95081e-05,5.46607e-08,-2.10555e-11,20323.9,17.1954], Tmin=(100,'K'), Tmax=(726.728,'K')), NASAPolynomial(coeffs=[4.91947,0.0223986,-1.16424e-05,2.34741e-09,-1.68432e-13,20063.2,6.2856], Tmin=(726.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
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
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,247.893,1569.59,1569.59,2285.13],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89822,0.00612407,0.000102608,-2.16319e-07,1.34206e-10,4932.65,9.22518], Tmin=(10,'K'), Tmax=(550.566,'K')), NASAPolynomial(coeffs=[3.92734,0.0284965,-1.98739e-05,6.49561e-09,-8.00611e-13,4587.15,5.99356], Tmin=(550.566,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(40.9583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCCDC1F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[C]=CC(F)=CF(893)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {7,S}
5 C u0 p0 c0 {3,S} {6,D} {8,S}
6 C u0 p1 c0 {5,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (97.7909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,194,682,905,1196,1383,3221,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.584831,'amu*angstrom^2'), symmetry=1, barrier=(13.4464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69408,0.0571226,-9.88748e-05,9.04918e-08,-3.21693e-11,11838.4,16.4177], Tmin=(100,'K'), Tmax=(822.3,'K')), NASAPolynomial(coeffs=[6.90619,0.0203029,-1.07944e-05,2.12507e-09,-1.48236e-13,11368.8,-5.35087], Tmin=(822.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.7909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsCsH) + group(CdJ2_singlet-Cds)"""),
)

species(
    label = 'F2(78)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[C]1C=C=C1(894)',
    structure = adjacencyList("""1 C u0 p0 c0 {3,S} {4,D} {5,S}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u0 p1 c0 {1,S} {2,S}
4 C u0 p0 c0 {1,D} {2,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (932.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,180,900.902,900.92,900.96,901.001,901.05,901.056,901.097],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0585,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.26862,0.0107893,1.39453e-06,-3.51803e-09,7.41548e-13,112103,16.0512], Tmin=(100,'K'), Tmax=(2156.13,'K')), NASAPolynomial(coeffs=[11.6575,0.00851458,-4.97673e-06,9.11216e-10,-5.71666e-14,106259,-31.4354], Tmin=(2156.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(932.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'CHF(40)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (138.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'FC1=C=C1(920)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (310.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0383,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04038,0.0168021,-5.40075e-06,-8.25794e-09,5.21779e-12,37381.9,7.19609], Tmin=(100,'K'), Tmax=(953.463,'K')), NASAPolynomial(coeffs=[8.6971,0.0035386,-1.00243e-06,1.81245e-10,-1.41139e-14,35827.4,-22.3199], Tmin=(953.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCddCF) + group(Cdd-CdsCds) + ring(Cyclopropadiene)"""),
)

species(
    label = 'F[C]1C[C]=C1F(899)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4 C u1 p0 c0 {1,S} {3,S} {5,S}
5 C u0 p0 c0 {2,S} {4,S} {6,D}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (195.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,271,519,563,612,1379,180,631.276,1553.96,1554.28,1554.66],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25213,0.0330749,-2.30923e-05,6.39674e-09,-2.96531e-13,23567.7,17.3148], Tmin=(100,'K'), Tmax=(1236.49,'K')), NASAPolynomial(coeffs=[10.8738,0.0117315,-5.14283e-06,1.00135e-09,-7.1479e-14,20935.1,-28.1336], Tmin=(1236.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(Cd-CdF1s)_ring) + radical(Cdj(Cs-CsHH)(Cd-CsF1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C]1C=[C]C1F(921)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (188.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,346,659,817,1284,2950,1000,180,931.927,1376.49,1376.77],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26921,0.0322639,-1.94956e-05,2.09376e-09,1.26536e-12,22712.1,17.5599], Tmin=(100,'K'), Tmax=(1148.6,'K')), NASAPolynomial(coeffs=[10.8624,0.0118522,-5.26379e-06,1.04481e-09,-7.59265e-14,20110.5,-27.8197], Tmin=(1148.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCdH)(F1s)(Cd-CdH)_ring) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1[C]=[C]C1F(922)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
5 C u1 p0 c0 {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (256.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([139,189,248,376,519,603,558,750,880,916,1147,1267,1279,1319,3102,3232,1830.96,1831.03],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48575,0.0257076,-5.85052e-07,-1.73047e-08,7.96289e-12,30916.4,16.221], Tmin=(100,'K'), Tmax=(1050.19,'K')), NASAPolynomial(coeffs=[10.5512,0.0119841,-5.25975e-06,1.07328e-09,-8.04988e-14,28285,-27.5476], Tmin=(1050.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1=C(F)[CH][CH]1(923)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (89.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,206,336,431,607,515,611,528,696,1312,1446,180,692.409,1476.48,1721.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3235.35,'J/mol'), sigma=(5.2451,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=505.35 K, Pc=50.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19216,0.0393646,-3.74e-05,1.78194e-08,-3.40083e-12,10888.1,15.3831], Tmin=(100,'K'), Tmax=(1250.71,'K')), NASAPolynomial(coeffs=[10.2921,0.0134598,-6.33194e-06,1.25922e-09,-9.07053e-14,8861.91,-25.5064], Tmin=(1250.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(Cd-CdF1s)(H)_ring) + radical(Csj(Cs-CdHH)(Cd-CdF1s)(H)_ring) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'FC1=[C]C=C1(901)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 C u0 p0 c0 {3,D} {4,S} {6,S}
3 C u0 p0 c0 {2,D} {5,S} {7,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (450.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,250,446,589,854,899,180,1020.27,1144.67,1154.18,1159.32,1654.21],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94113,0.00329487,7.45308e-05,-1.35534e-07,7.37088e-11,54134.3,9.75919], Tmin=(10,'K'), Tmax=(606.155,'K')), NASAPolynomial(coeffs=[2.49265,0.0270597,-1.94333e-05,6.4748e-09,-8.0767e-13,54048.9,13.8693], Tmin=(606.155,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(450.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1D[C]CDC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1[C]=C=C1(924)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3 C u0 p0 c0 {2,S} {5,D} {7,S}
4 C u1 p0 c0 {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (542.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,180,180,180,1025.5,1645.73,1654.56],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87823,0.0310553,-5.3966e-05,6.09904e-08,-2.60864e-11,65277.8,12.6061], Tmin=(100,'K'), Tmax=(801.591,'K')), NASAPolynomial(coeffs=[0.805769,0.0254758,-1.37323e-05,2.75077e-09,-1.94872e-13,66121.5,25.3368], Tmin=(801.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cds_S)"""),
)

species(
    label = 'FC1=[C]C=C1F(902)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (271.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2950,1000,250,446,589,854,899,678.518,1615.58,2023.1],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32993,0.0305922,-1.85647e-05,-1.12734e-09,3.23424e-12,32712.4,15.7183], Tmin=(100,'K'), Tmax=(1029.61,'K')), NASAPolynomial(coeffs=[11.4765,0.00752368,-3.11766e-06,6.29656e-10,-4.75468e-14,30168.2,-31.8841], Tmin=(1029.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cdj(Cd-CdF1s)(Cd-CdH)_ring) + longDistanceInteraction_cyclic(Cds(F)-Cds(F)) + longDistanceInteraction_cyclic(Cds(F)-Cds(F))"""),
)

species(
    label = 'FC1=C=[C]C1F(925)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4 C u0 p0 c0 {2,S} {3,S} {6,D}
5 C u1 p0 c0 {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (366.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,145,326,398,834,1303,180,875.208,2110.4],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.53676,0.019415,-9.45815e-06,1.34419e-09,-2.85511e-14,44044.1,8.19212], Tmin=(100,'K'), Tmax=(2912.91,'K')), NASAPolynomial(coeffs=[39.4018,-0.0183175,4.74849e-06,-7.11702e-10,4.52891e-14,19428.7,-204.675], Tmin=(2912.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CdCddCF) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Cds_S) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'C1=C=CC=1(903)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,D} {6,S}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (652.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,180,202.072,1081.32,1081.4,1081.42,1081.48,1081.54,1082.1],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0585,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.23155,0.0114051,1.27491e-05,-2.36761e-08,9.44414e-12,78559.5,7.91396], Tmin=(100,'K'), Tmax=(996.209,'K')), NASAPolynomial(coeffs=[7.20433,0.00882819,-3.50925e-06,6.80706e-10,-5.02341e-14,77104.3,-14.5683], Tmin=(996.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
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
    E0 = (32.8384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (41.1227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (257.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (255.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (284.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (509.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (298.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (195.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (224.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (187.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (187.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (169.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (10.1932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (56.4353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (320.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (263.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (63.0693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (138.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (188.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (273.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-4.04871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (137.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (176.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (410.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (318.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (205.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (45.6362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (283.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (226.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (28.4063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (73.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-13.2987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (731.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (413.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (322.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (451.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (203.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (31.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (24.3358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (128.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-38.4125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (331.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (423.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (291.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (386.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (459.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (231.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['C2HF(58)', 'C2HF(58)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;Cdsingle_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C2HF(58)', '[CH]=[C]F(252)'],
    products = ['[CH]=C(F)C=[C]F(246)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0243,'m^3/(mol*s)'), n=2.41, Ea=(3.02832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H_Ext-3R-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(5)', '[CH]=C(F)C#CF(308)'],
    products = ['[CH]=C(F)C=[C]F(246)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(763,'m^3/(mol*s)'), n=1.64, Ea=(7.62199,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', 'C#CC=[C]F(309)'],
    products = ['[CH]=C(F)C=[C]F(246)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.2366,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]F(252)', '[CH]=[C]F(252)'],
    products = ['[CH]=C(F)C=[C]F(246)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -36.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['HF(38)', '[CH]=C=C=[C]F(310)'],
    products = ['[CH]=C(F)C=[C]F(246)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(262.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['[CH]C(F)=C=CF(311)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.50921e+08,'s^-1'), n=1.36125, Ea=(162.666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['[CH2]C(F)=C=[C]F(293)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['[CH]=C=C[C](F)F(282)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(155.088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(F)C=[C]F(246)'],
    products = ['F[C]C=C=CF(306)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(154.87,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['FC1(F)C=C=C1(572)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(284.614,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['FC1=[C]C(F)[CH]1(573)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]1[C]=C(F)C1(574)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'FC1=C[C]=C1(575)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'FC1=[C]C(F)=C1(576)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.03232e+06,'m^3/(mol*s)'), n=-0.00929868, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Sp-6R!H-5R!H',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Sp-6R!H-5R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.03468e+20,'s^-1'), n=-2.18552, Ea=(85.6267,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5822770821690705, var=6.058307656543806, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['FC1=C[C]C1F(578)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(151.739,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', 'FC1=C=C=C1(579)'],
    products = ['FC1=CC(F)=C1(307)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(191.534,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1(F)[C]C=C1(862)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(296.27,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]=C(F)C=CF(863)'],
    products = ['FC1=CC(F)[C]1(577)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1=CC=C1F(304)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(160.15,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FC1=CC(F)[C]1(577)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(199.312,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'FC1=C=C[C]1(865)'],
    products = ['FC1=CC(F)[C]1(577)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(142.189,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'FC1[C]C#C1(866)'],
    products = ['FC1=CC(F)[C]1(577)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(386.715,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]=C=CC(F)F(874)'],
    products = ['[C]=C(F)C=CF(863)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.45932e+11,'s^-1'), n=0.63878, Ea=(287.201,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]=C(F)C=CF(863)'],
    products = ['FC#CC=CF(875)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(172.098,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', '[C]=C(F)C#C(876)'],
    products = ['[C]=C(F)C=CF(863)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(322.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[C]=C=C=CF(877)'],
    products = ['[C]=C(F)C=CF(863)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(191.67,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[C]=CC(F)=CF(893)'],
    products = ['FC1=C[C]C1F(578)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FC1=C[C]C1F(578)'],
    products = ['FC1=CC=C1F(304)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.03468e+20,'s^-1'), n=-2.18552, Ea=(86.5356,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5822770821690705, var=6.058307656543806, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FC1=C[C]C1F(578)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.26313e-09,'s^-1'), n=5.53373, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.525240051516825, var=52.082221049180056, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_N-Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_N-Sp-5R!H-1C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F2(78)', '[C]1C=C=C1(894)'],
    products = ['FC1=C[C]C1F(578)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.000237309,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'FC1=C=C[C]1(865)'],
    products = ['FC1=C[C]C1F(578)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(145.285,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'FC1[C]C#C1(866)'],
    products = ['FC1=C[C]C1F(578)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(391.413,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CHF(40)', 'FC1=C=C1(920)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(193.701,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC1=C=CC1F(864)'],
    products = ['FC1=CC=C1F(304)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(266.096,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction36',
    reactants = ['F[C]1C[C]=C1F(899)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['F[C]1C=[C]C1F(921)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['FC1[C]=[C]C1F(922)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['FC1=C(F)[CH][CH]1(923)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['F(37)', 'FC1=[C]C=C1(901)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction41',
    reactants = ['F(37)', 'FC1[C]=C=C1(924)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(5)', 'FC1=[C]C=C1F(902)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(9.3466e+06,'m^3/(mol*s)'), n=0.408973, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_N-4R!H-inRing',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_N-4R!H-inRing"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(5)', 'FC1=C=[C]C1F(925)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.37186e+06,'m^3/(mol*s)'), n=0.39882, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction44',
    reactants = ['F2(78)', 'C1=C=CC=1(903)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.000237309,'m^3/(mol*s)'), n=2.63647, Ea=(7.3193,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction45',
    reactants = ['HF(38)', 'FC1=C=C=C1(579)'],
    products = ['FC1=C=CC1F(864)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5353.26,'m^3/(mol*s)'), n=0.732206, Ea=(235.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #46',
    isomers = [
        '[CH]=C(F)C=[C]F(246)',
        'FC1=CC(F)=C1(307)',
        'FC1=CC(F)[C]1(577)',
        '[C]=C(F)C=CF(863)',
        'FC1=C[C]C1F(578)',
        'FC1=C=CC1F(864)',
    ],
    reactants = [
        ('C2HF(58)', 'C2HF(58)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #46',
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

