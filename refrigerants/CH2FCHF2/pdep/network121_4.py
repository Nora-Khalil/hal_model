species(
    label = '[O]OC1(F)[CH]C(F)=C1(581)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-149.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,268.418,268.442,268.452,268.743,268.796],'cm^-1')),
        HinderedRotor(inertia=(0.910361,'amu*angstrom^2'), symmetry=1, barrier=(46.8537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3903.93,'J/mol'), sigma=(6.3222,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=609.78 K, Pc=35.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02254,0.0615362,-6.13767e-05,2.96625e-08,-5.63524e-12,-17846.6,19.5807], Tmin=(100,'K'), Tmax=(1274.42,'K')), NASAPolynomial(coeffs=[15.6097,0.0157519,-7.4884e-06,1.47287e-09,-1.05362e-13,-21564.6,-54.3311], Tmin=(1274.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'FC1=CC2(F)OOC12(939)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-248.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,323,467,575,827,1418,556.295,556.354,556.401,556.506,556.667,556.693,556.734,556.988,556.993],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3642.96,'J/mol'), sigma=(5.89185,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.02 K, Pc=40.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37505,0.0476133,-2.48125e-05,-2.36513e-09,3.49961e-12,-29797.1,16.5191], Tmin=(100,'K'), Tmax=(1189.18,'K')), NASAPolynomial(coeffs=[15.9615,0.016946,-9.33452e-06,1.96667e-09,-1.45909e-13,-34567,-61.8485], Tmin=(1189.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = '[O]C1C(F)=CC1([O])F(964)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-201.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,323,467,575,827,1418,561.027,561.027,561.027,561.027,561.027,561.027,561.027,561.027,561.027],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4271.9,'J/mol'), sigma=(6.7089,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=667.26 K, Pc=32.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14284,0.052484,-3.48088e-05,1.04848e-09,4.54051e-12,-24160.8,22.6183], Tmin=(100,'K'), Tmax=(1035.18,'K')), NASAPolynomial(coeffs=[17.0604,0.0114786,-5.098e-06,1.04607e-09,-7.9314e-14,-28554.8,-60.0307], Tmin=(1035.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(F)C([O])C(=O)F(1465)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {4,D} {5,S}
8  C u1 p0 c0 {6,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-230.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,246,474,533,1155,486,617,768,1157,1926,3120,650,792.5,1650,294.366,294.653,4000],'cm^-1')),
        HinderedRotor(inertia=(0.164285,'amu*angstrom^2'), symmetry=1, barrier=(10.129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0169845,'amu*angstrom^2'), symmetry=1, barrier=(26.5656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43475,0.0597448,-7.34796e-05,4.72429e-08,-1.22353e-11,-27673.5,27.902], Tmin=(100,'K'), Tmax=(935.911,'K')), NASAPolynomial(coeffs=[10.7146,0.0200838,-9.9145e-06,1.96446e-09,-1.40571e-13,-29410.5,-16.2532], Tmin=(935.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CdCsCdF) + group(COCsFO) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.0012157,5.31615e-06,-4.8944e-09,1.45844e-12,-1038.59,4.68369], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15383,0.00167803,-7.69968e-07,1.51274e-10,-1.08781e-14,-1040.82,6.16752], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'O=O(207)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 O u0 p2 c0 {1,D}
"""),
    E0 = (85.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1857.18,'J/mol'), sigma=(4.34667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=290.09 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.0012157,5.31615e-06,-4.8944e-09,1.45844e-12,10304.5,4.68369], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15383,0.00167803,-7.69968e-07,1.51274e-10,-1.08781e-14,10302.3,6.16752], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]O[C](F)C1C=C1F(936)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,D}
8  C u1 p0 c0 {2,S} {3,S} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (94.4906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,323,467,575,827,1418,395,473,707,1436,180,180,180,1379.2,1379.2,1379.24],'cm^-1')),
        HinderedRotor(inertia=(0.177033,'amu*angstrom^2'), symmetry=1, barrier=(4.07034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177161,'amu*angstrom^2'), symmetry=1, barrier=(4.07328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22454,0.0665684,-0.000101708,8.68776e-08,-2.99134e-11,11459.2,25.6449], Tmin=(100,'K'), Tmax=(781.525,'K')), NASAPolynomial(coeffs=[8.21168,0.0255525,-1.29006e-05,2.5194e-09,-1.76388e-13,10527.6,-5.31502], Tmin=(781.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.4906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd-Cd-Cs(C-F)) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC1=CC(F)(F)[CH]1(937)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {7,S} {8,D}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u0 p0 c0 {5,S} {6,D} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-29.6977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,248,333,466,604,684,796,1061,1199,2750,3150,900,1100,180,180,896.311,1420.31,1420.42,1420.43,1420.44,4000],'cm^-1')),
        HinderedRotor(inertia=(0.30955,'amu*angstrom^2'), symmetry=1, barrier=(7.11716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60742,0.0476672,-4.13733e-05,1.76967e-08,-3.00191e-12,-3481.64,25.1617], Tmin=(100,'K'), Tmax=(1411.68,'K')), NASAPolynomial(coeffs=[13.1887,0.0148516,-6.50455e-06,1.22991e-09,-8.57382e-14,-6751.44,-34.7042], Tmin=(1411.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.6977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCCFF) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(ROOJ) + radical(CCJCO)"""),
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
    label = '[O]C1(F)[CH]C(F)=C1(938)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u0 p0 c0 {4,S} {7,D} {9,S}
7 C u0 p0 c0 {2,S} {5,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-142.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,271,519,563,612,1379,341.278,342.214,342.774,343.463,344.342,344.65],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46585,0.045777,-2.52862e-05,-5.30751e-09,6.19495e-12,-17033.1,16.276], Tmin=(100,'K'), Tmax=(1028.49,'K')), NASAPolynomial(coeffs=[15.6289,0.011289,-5.024e-06,1.02852e-09,-7.78941e-14,-21035.7,-57.7457], Tmin=(1028.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'FC1=CC2(F)OC12(940)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 C u0 p0 c0 {4,S} {7,D} {9,S}
7 C u0 p0 c0 {2,S} {5,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-250.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76629,0.0397829,-1.86271e-05,-6.18262e-09,4.91723e-12,-30036.3,14.5526], Tmin=(100,'K'), Tmax=(1112.21,'K')), NASAPolynomial(coeffs=[14.3568,0.012736,-6.74177e-06,1.43388e-09,-1.08173e-13,-33964.7,-52.5985], Tmin=(1112.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_3_4_ene_1)"""),
)

species(
    label = 'HO2(11)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.49012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[O]OC1(F)C2[C](F)C21(941)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-19.0462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67584,0.0561009,-5.96508e-05,3.3299e-08,-7.7483e-12,-2211.15,19.3111], Tmin=(100,'K'), Tmax=(1013.28,'K')), NASAPolynomial(coeffs=[9.63708,0.024673,-1.31264e-05,2.68896e-09,-1.9603e-13,-3824.54,-19.2022], Tmin=(1013.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.0462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCCFO) + group(CsCsCsFH) + polycyclic(s2_3_3_ane) + radical(ROOJ) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FC12[CH]C(F)([CH]1)OO2(942)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-61.3995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52271,0.0425522,-8.32105e-06,-1.9032e-08,8.96997e-12,-7285.51,17.5608], Tmin=(100,'K'), Tmax=(1131.21,'K')), NASAPolynomial(coeffs=[15.684,0.0182186,-1.01874e-05,2.18378e-09,-1.64439e-13,-12136.3,-59.7845], Tmin=(1131.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.3995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + polycyclic(s3_4_5_ane) + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = 'F[C]1[CH]C2(F)OOC12(943)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {5,S} {8,S} {10,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (24.0771,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,3150,900,1100,212,367,445,1450,792.942,792.942,792.942,792.942,792.943,792.943,792.943,792.943,792.943,792.943],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01799,0.0472587,-3.27286e-05,9.66921e-09,-1.11055e-12,2962.68,18.0674], Tmin=(100,'K'), Tmax=(1964.1,'K')), NASAPolynomial(coeffs=[16.3603,0.0180501,-1.04219e-05,2.09781e-09,-1.46839e-13,-2671.3,-60.8074], Tmin=(1964.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.0771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + polycyclic(s2_4_4_ane) + radical(CCJCOOH) + radical(CsCsCsF1s)"""),
)

species(
    label = '[CH]=C(F)C=C(F)O[O](944)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {6,D} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {3,S} {5,D}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (44.4727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,250,446,589,854,899,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.393707,'amu*angstrom^2'), symmetry=1, barrier=(9.0521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393344,'amu*angstrom^2'), symmetry=1, barrier=(9.04375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889666,0.0760107,-0.000129009,1.14687e-07,-3.96583e-11,5453.51,26.1014], Tmin=(100,'K'), Tmax=(828.651,'K')), NASAPolynomial(coeffs=[8.87997,0.0246153,-1.2758e-05,2.48311e-09,-1.72016e-13,4569.61,-8.28846], Tmin=(828.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.4727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(CdCFO) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cdj(Cd-CdF1s)(H))"""),
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
    label = '[O]OC1=CC(F)=C1(945)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {2,S} {5,S} {6,D}
5 C u0 p0 c0 {4,S} {7,D} {8,S}
6 C u0 p0 c0 {4,D} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (277.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,280,518,736,852,873,248.088,248.208,248.638,249.113,1919.05,1919.96,1920.46,3443],'cm^-1')),
        HinderedRotor(inertia=(0.287543,'amu*angstrom^2'), symmetry=1, barrier=(12.6568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72327,0.0531044,-7.25721e-05,5.35399e-08,-1.58478e-11,33454.4,20.3365], Tmin=(100,'K'), Tmax=(825.594,'K')), NASAPolynomial(coeffs=[9.0366,0.0176726,-8.19922e-06,1.56063e-09,-1.08397e-13,32246.8,-13.5447], Tmin=(825.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + ring(Cd-Cd-Cd-Cd(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1(F)C=C=C1(946)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {7,D} {8,S}
6 C u0 p0 c0 {4,S} {7,D} {9,S}
7 C u0 p0 c0 {5,D} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (310.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,3150,900,1100,180,180,180,180,949.25,1646.19,1651.09],'cm^-1')),
        HinderedRotor(inertia=(1.66094,'amu*angstrom^2'), symmetry=1, barrier=(38.1884,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61581,0.0601179,-0.000102199,9.82056e-08,-3.70641e-11,37460.4,17.8954], Tmin=(100,'K'), Tmax=(799.523,'K')), NASAPolynomial(coeffs=[5.07374,0.0283104,-1.53064e-05,3.05656e-09,-2.15843e-13,37371.1,4.88633], Tmin=(799.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(ROOJ)"""),
)

species(
    label = 'F[C]1[CH]C(F)=C1(615)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (86.8055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,271,519,563,612,1379,180,180,1104.81,1106.27,1801.34],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06134,0.0383117,-3.38918e-05,1.45559e-08,-2.45109e-12,10513.6,15.9062], Tmin=(100,'K'), Tmax=(1431.45,'K')), NASAPolynomial(coeffs=[12.0907,0.0102863,-4.52462e-06,8.78888e-10,-6.24417e-14,7642.3,-36.0769], Tmin=(1431.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.8055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(Cd-CdH)_ring) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring)"""),
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
    label = '[O]OC1=C=C(F)[CH]1(947)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {2,S} {5,S} {7,D}
5 C u1 p0 c0 {4,S} {6,S} {8,S}
6 C u0 p0 c0 {1,S} {5,S} {7,D}
7 C u0 p0 c0 {4,D} {6,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (495.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,1000,132,421,945,1169,1268,180,180,180,180,1244.42,1245.93,1248.57],'cm^-1')),
        HinderedRotor(inertia=(0.203759,'amu*angstrom^2'), symmetry=1, barrier=(4.68483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75137,0.0548661,-8.90412e-05,7.9964e-08,-2.86758e-11,59617.6,20.2984], Tmin=(100,'K'), Tmax=(783.379,'K')), NASAPolynomial(coeffs=[6.87962,0.0215188,-1.1475e-05,2.28365e-09,-1.61248e-13,59033.9,-1.78781], Tmin=(783.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(CdCddCF) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC1(F)[C]=C=C1(948)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {7,D} {8,S}
6 C u1 p0 c0 {4,S} {7,D}
7 C u0 p0 c0 {5,D} {6,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (548.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,1000,180,180,180,180,848.035,851.095],'cm^-1')),
        HinderedRotor(inertia=(0.149384,'amu*angstrom^2'), symmetry=1, barrier=(3.43464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5568,0.0650164,-0.000127944,1.28457e-07,-4.82537e-11,66064.8,19.2741], Tmin=(100,'K'), Tmax=(840.519,'K')), NASAPolynomial(coeffs=[4.08905,0.0274515,-1.53725e-05,3.05478e-09,-2.12908e-13,66540.4,12.8586], Tmin=(840.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC1(F)[C]=C(F)C1(949)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-13.7343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,3150,900,1100,246,474,533,1155,465.805,466.62,468.299,468.408,468.508,469.41],'cm^-1')),
        HinderedRotor(inertia=(0.000772044,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46796,0.0590116,-6.69781e-05,3.99845e-08,-9.73019e-12,-1563.41,21.2655], Tmin=(100,'K'), Tmax=(986.864,'K')), NASAPolynomial(coeffs=[10.5465,0.0222134,-1.10452e-05,2.19908e-09,-1.57945e-13,-3355.25,-22.4132], Tmin=(986.864,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.7343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'OOC1(F)[C]=C(F)[CH]1(950)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u1 p0 c0 {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-48.8223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.941702,0.0671657,-7.24622e-05,3.7943e-08,-7.86338e-12,-5761.73,19.5874], Tmin=(100,'K'), Tmax=(1165.68,'K')), NASAPolynomial(coeffs=[15.201,0.0182359,-9.50002e-06,1.93459e-09,-1.40864e-13,-9086.13,-51.3917], Tmin=(1165.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.8223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(C=CCJCO) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'FOO[C]1[CH]C(F)=C1(951)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {2,S} {3,S}
5  C u1 p0 c0 {3,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (229.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,277,555,632,2750,3150,900,1100,271,519,563,612,1379,526.961,526.963,526.967,526.968,526.968],'cm^-1')),
        HinderedRotor(inertia=(0.242903,'amu*angstrom^2'), symmetry=1, barrier=(47.8664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24291,'amu*angstrom^2'), symmetry=1, barrier=(47.8666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41157,0.055115,-4.7738e-05,1.95688e-08,-3.20082e-12,27752.3,20.6038], Tmin=(100,'K'), Tmax=(1438.08,'K')), NASAPolynomial(coeffs=[14.4963,0.0187198,-9.77551e-06,1.96997e-09,-1.4137e-13,23988.9,-47.2762], Tmin=(1438.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd(F)-Cd-Cs-Cs) + radical(C2CsJOO) + radical(C=CCJCO)"""),
)

species(
    label = '[O]O[C]1C=C(F)C1F(952)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {3,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-19.6738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,174,267,591,721,1107,1278,1348,3273,323,467,575,827,1418,2950,1000,180,1358.58,1359.67,1361.92,1363.82],'cm^-1')),
        HinderedRotor(inertia=(0.333706,'amu*angstrom^2'), symmetry=1, barrier=(7.67256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93229,0.0472551,-4.08128e-05,1.81014e-08,-3.3382e-12,-2293.41,23.6488], Tmin=(100,'K'), Tmax=(1253.43,'K')), NASAPolynomial(coeffs=[9.79843,0.0221522,-1.07717e-05,2.12328e-09,-1.51324e-13,-4265.34,-16.0776], Tmin=(1253.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(C2CsJOOH) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'FOOC1(F)[CH][C]=C1(953)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {2,S} {3,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (254.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,277,555,632,1380,1390,370,380,2900,435,2750,3150,900,1100,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.93967,'amu*angstrom^2'), symmetry=1, barrier=(44.5968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94074,'amu*angstrom^2'), symmetry=1, barrier=(44.6214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01737,0.0633608,-6.29435e-05,2.995e-08,-5.63057e-12,30775.8,19.9049], Tmin=(100,'K'), Tmax=(1278.9,'K')), NASAPolynomial(coeffs=[15.7265,0.0173553,-8.98441e-06,1.82222e-09,-1.32138e-13,27013.5,-54.6767], Tmin=(1278.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sFO) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(C=CCJCO) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[O]OC1(F)C=[C]C1F(954)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (20.6394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,2950,1000,180,180,752.336,752.395],'cm^-1')),
        HinderedRotor(inertia=(0.108928,'amu*angstrom^2'), symmetry=1, barrier=(43.7635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57544,0.0551226,-5.51713e-05,2.8131e-08,-5.84039e-12,2568.15,21.123], Tmin=(100,'K'), Tmax=(1143.82,'K')), NASAPolynomial(coeffs=[11.3138,0.0210674,-1.0512e-05,2.10186e-09,-1.51366e-13,340.347,-27.1675], Tmin=(1143.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.6394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCCFO) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1=CC2(OO2)C1F(961)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-221.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([284,328,853,1146,1135,1297,3239,2950,1000,323,467,575,827,1418,514.234,514.848,515.067,515.133,515.347,515.535,515.846,515.874,516.035,516.689],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.952518,0.0513104,-1.49392e-05,-3.07211e-08,1.85052e-11,-26554.3,17.9795], Tmin=(100,'K'), Tmax=(964.911,'K')), NASAPolynomial(coeffs=[21.2852,0.0047466,-1.19788e-06,3.02741e-10,-3.061e-14,-32234.3,-88.4876], Tmin=(964.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + polycyclic(s1_3_4_ene)"""),
)

species(
    label = 'FC1(F)C=C2OOC21(962)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-197.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,274,345,380,539,705,1166,1213,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6909,0.0426832,-2.2114e-05,-3.93846e-10,2.07861e-12,-23616.1,18.4544], Tmin=(100,'K'), Tmax=(1273.93,'K')), NASAPolynomial(coeffs=[14.3701,0.0174185,-9.49376e-06,1.96511e-09,-1.43292e-13,-28027,-50.4183], Tmin=(1273.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFF) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(12dioxetane)"""),
)

species(
    label = 'F[C]1C=C(F)[CH]OO1(963)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {6,D} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,D} {7,S}
7  C u1 p0 c0 {4,S} {6,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-149.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,271,519,563,612,1379,3025,407.5,1350,352.5,686.578,686.578,686.578,686.578,686.578,686.578,686.578,686.578,686.578,686.578,686.579],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68693,0.0312264,3.57095e-05,-7.72132e-08,3.32851e-11,-17927.7,18.6024], Tmin=(100,'K'), Tmax=(970.013,'K')), NASAPolynomial(coeffs=[19.2034,0.00748418,-2.55937e-06,6.22117e-10,-5.71252e-14,-23607.2,-77.1301], Tmin=(970.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFHO) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(36dihydro12dioxin) + radical(Cs_P) + radical(C=CCJO)"""),
)

species(
    label = '[O]OC1[C](F)C=C1F(580)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {2,S} {5,S} {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-83.7573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,180,180,293.128,1359.42,1359.74,1359.79,1359.88],'cm^-1')),
        HinderedRotor(inertia=(1.14334,'amu*angstrom^2'), symmetry=1, barrier=(26.2877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26433,0.0521615,-4.29954e-05,1.47154e-08,-1.26522e-12,-9968.42,24.5603], Tmin=(100,'K'), Tmax=(1133.33,'K')), NASAPolynomial(coeffs=[15.1191,0.0137697,-6.0895e-06,1.18641e-09,-8.53912e-14,-13783.6,-46.992], Tmin=(1133.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.7573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(ROOJ) + radical(CsCdCsF1s)"""),
)

species(
    label = 'F[C]=CC1(F)[CH]OO1(965)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u1 p0 c0 {4,S} {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {2,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (55.5743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,464.769,464.985,465.058,465.255,465.266,465.283,465.529],'cm^-1')),
        HinderedRotor(inertia=(0.000778836,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4877,0.0603484,-7.81999e-05,5.70322e-08,-1.72581e-11,6769.95,22.8647], Tmin=(100,'K'), Tmax=(796.761,'K')), NASAPolynomial(coeffs=[8.36637,0.0258145,-1.31843e-05,2.63115e-09,-1.88296e-13,5673.84,-8.75797], Tmin=(796.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.5743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCsJOO) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=C(F)C1OO[C]1F(966)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {1,S} {4,S} {5,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (73.0536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,246,474,533,1155,3120,650,792.5,1650,180,180,180,1035.36,1035.4,1035.4,1035.41,1035.41,3737.04],'cm^-1')),
        HinderedRotor(inertia=(0.157489,'amu*angstrom^2'), symmetry=1, barrier=(3.62099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61723,0.0542419,-5.67991e-05,3.05369e-08,-6.65842e-12,8870.62,25.6983], Tmin=(100,'K'), Tmax=(1096.5,'K')), NASAPolynomial(coeffs=[11.0395,0.0198694,-9.77762e-06,1.94783e-09,-1.40133e-13,6804.33,-20.6268], Tmin=(1096.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.0536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = 'FC1=COOC(F)=C1(967)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  C u0 p0 c0 {4,S} {6,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-254.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,280,518,736,852,873,326,540,652,719,1357,180,180,1192.72,1192.78,1193.06,1193.08,1193.14,1193.25,1193.37,1193.47],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10305,0.0406104,-2.38555e-05,5.52708e-09,-4.23174e-13,-30543.9,18.0078], Tmin=(100,'K'), Tmax=(1927.43,'K')), NASAPolynomial(coeffs=[17.7723,0.0138948,-7.58039e-06,1.4598e-09,-9.8227e-14,-37662,-70.6656], Tmin=(1927.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(1,3-Cyclohexadiene)"""),
)

species(
    label = 'FC1[CH]C2(F)OO[C]12(968)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (20.2264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,259,529,569,1128,1321,1390,3140,2950,1000,694.38,694.677,694.84,694.859,695.415,695.524,695.564,695.672,695.699],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14618,0.0459413,-3.00887e-05,8.14743e-09,-8.40894e-13,2493.5,15.7014], Tmin=(100,'K'), Tmax=(2209.68,'K')), NASAPolynomial(coeffs=[19.2066,0.0150586,-9.12473e-06,1.82262e-09,-1.25323e-13,-5046.2,-80.1319], Tmin=(2209.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CCJCOOH)"""),
)

species(
    label = 'F[C]1CC2(F)OO[C]12(969)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {8,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (10.2856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([316,385,515,654,689,1295,2750,3150,900,1100,212,367,445,1450,830.071,830.318,830.397,830.41,830.422,830.564,830.612,830.687,830.714,830.718],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64904,0.0427419,-2.60229e-05,6.22733e-09,-5.51256e-13,1271.99,14.2565], Tmin=(100,'K'), Tmax=(2659.87,'K')), NASAPolynomial(coeffs=[25.9948,0.00763381,-6.22418e-06,1.26499e-09,-8.48486e-14,-11147.3,-121.212], Tmin=(2659.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.2856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + polycyclic(s2_4_4_ane) + radical(C2CsJOO) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FC1=C[C]2OOC21(970)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u0 p0 c0 {5,S} {6,D} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (182.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,323,467,575,827,1418,819.899,819.915,819.917,819.924,819.924,819.925,819.926,819.928,819.931,819.94,819.95,819.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24337,0.0284262,8.61439e-06,-2.86093e-08,1.12652e-11,22041.7,17.3747], Tmin=(100,'K'), Tmax=(1102.8,'K')), NASAPolynomial(coeffs=[11.9075,0.0174398,-9.17679e-06,1.93459e-09,-1.44964e-13,18446.7,-36.8297], Tmin=(1102.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(CdCsCdF) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOO)"""),
)

species(
    label = 'FC12C=[C]C1OO2(971)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 C u0 p0 c0 {4,S} {7,D} {9,S}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (212.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,600.775,602.021,602.064,602.889,603.065,603.084,603.517,603.719,603.852,604.829,606.742],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90063,0.0373054,-1.2475e-05,-9.21548e-09,4.97158e-12,25626.6,15.0075], Tmin=(100,'K'), Tmax=(1182.18,'K')), NASAPolynomial(coeffs=[13.4827,0.0163611,-9.04977e-06,1.90788e-09,-1.41524e-13,21613.2,-48.1999], Tmin=(1182.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutene-vinyl)"""),
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
    label = 'FC1=CC2(F)OO[C]12(972)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u1 p0 c0 {4,S} {5,S} {8,S}
7 C u0 p0 c0 {5,S} {8,D} {9,S}
8 C u0 p0 c0 {2,S} {6,S} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-61.9787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,271,519,563,612,1379,563.979,563.98,563.992,564.069,564.133,564.247,564.388,564.602],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01417,0.045025,-3.24803e-05,9.94643e-09,-1.17886e-12,-7384.98,15.0075], Tmin=(100,'K'), Tmax=(1932.46,'K')), NASAPolynomial(coeffs=[16.6506,0.0147295,-8.96495e-06,1.83416e-09,-1.29406e-13,-13041.9,-65.2476], Tmin=(1932.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.9787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + polycyclic(s2_4_4_ene_1) + radical(C2CsJOO)"""),
)

species(
    label = 'FC1=[C]C2(F)OOC12(973)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6 C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7 C u0 p0 c0 {2,S} {5,S} {8,D}
8 C u1 p0 c0 {6,S} {7,D}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (3.88367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,246,474,533,1155,549.196,549.199,549.199,549.201,549.202,549.203,549.204,549.208,549.212],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58137,0.0494147,-3.98254e-05,1.41196e-08,-1.94568e-12,556.87,16.2496], Tmin=(100,'K'), Tmax=(1708.1,'K')), NASAPolynomial(coeffs=[17.1503,0.0129556,-7.80814e-06,1.62335e-09,-1.16704e-13,-4761.79,-67.1969], Tmin=(1708.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.88367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CdCsCdF) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'FC1[C]C2(F)OOC12(974)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
8  C u0 p1 c0 {6,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (75.1403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,1380,1390,370,380,2900,435,353,444,1253,3145,503.275,503.278,503.28,503.28,503.28,503.281,503.281,503.281,503.281,503.282,503.282,503.293],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54239,0.0513201,-3.51889e-05,1.0202e-08,-1.13059e-12,9127.73,24.5604], Tmin=(100,'K'), Tmax=(2107.59,'K')), NASAPolynomial(coeffs=[21.6856,0.013091,-7.98122e-06,1.59587e-09,-1.09763e-13,636.845,-87.637], Tmin=(2107.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.1403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFH) + group(CsJ2_singlet-CsH) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'FC1[C]C2OOC12F(975)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u0 p2 c0 {3,S} {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
8  C u0 p1 c0 {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (100.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([333,384,448,608,1254,1480,2950,1000,353,444,1253,3145,600.7,600.702,600.702,600.713,600.719,600.719,600.726,600.728,600.729,600.731,600.735,600.736],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41929,0.0512318,-3.51613e-05,1.03328e-08,-1.16199e-12,12147.5,26.3235], Tmin=(100,'K'), Tmax=(2091.91,'K')), NASAPolynomial(coeffs=[21.814,0.0122349,-7.19893e-06,1.42155e-09,-9.70394e-14,3614.63,-87.1221], Tmin=(2091.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCCFH) + group(CsJ2_singlet-CsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'FC1=CC2=C1OO2(976)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {5,D} {6,S}
5 C u0 p0 c0 {3,S} {4,D} {7,S}
6 C u0 p0 c0 {4,S} {7,D} {8,S}
7 C u0 p0 c0 {1,S} {5,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (435.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,280,518,736,852,873,180,180,180,256.366,1388.64,1388.88,1389.76,1389.87,1389.94,1391.84,1393.66],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29868,0.0442337,-7.25845e-05,7.87401e-08,-3.41758e-11,52472.7,15.6432], Tmin=(100,'K'), Tmax=(749.133,'K')), NASAPolynomial(coeffs=[1.7902,0.032016,-1.82209e-05,3.75231e-09,-2.71132e-13,52967.9,20.7461], Tmin=(749.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + polycyclic(s2_4_4_diene_1_m)"""),
)

species(
    label = 'FC12C=C=C1OO2(977)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,D}
6 C u0 p0 c0 {4,S} {7,D} {8,S}
7 C u0 p0 c0 {5,D} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (346.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,689.329,689.33,689.33,689.33,689.33,689.331,689.331,689.332,689.332,689.333],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94025,0.0284906,-1.53689e-05,2.53487e-09,-8.92946e-14,41603.7,7.10815], Tmin=(100,'K'), Tmax=(2602.4,'K')), NASAPolynomial(coeffs=[34.4545,-0.00798926,-3.4926e-07,2.26052e-10,-1.53254e-14,22192.6,-176.071], Tmin=(2602.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsCCFO) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(cyclobutadiene_13) + ring(12dioxetane)"""),
)

species(
    label = 'FC1=C=C2OOC12(978)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5 C u0 p0 c0 {3,S} {4,S} {7,D}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u0 p0 c0 {5,D} {6,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (397.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,145,326,398,834,1303,180,180,296.008,1292.42,1293.17,1294.2,1295.14,1297.28,1299.35,1299.46,1300.41],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39646,0.0276662,-1.41735e-05,2.00555e-09,1.58983e-14,47778.8,13.1925], Tmin=(100,'K'), Tmax=(2270.01,'K')), NASAPolynomial(coeffs=[22.2229,0.00368778,-4.40527e-06,9.21346e-10,-6.12331e-14,36862.3,-98.2863], Tmin=(2270.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(CdCddCF) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(cyclobutadiene_13) + ring(12dioxetane)"""),
)

species(
    label = 'FC12C#CC1OO2(979)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 C u0 p0 c0 {4,S} {7,T}
7 C u0 p0 c0 {5,S} {6,T}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (7.16824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,397.803,398.93,400.735,401.109,401.204,401.221,402.458,402.51,402.556,403.269],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56298,0.0463651,-3.88567e-05,1.42418e-08,-1.69472e-12,955.935,-2.88482], Tmin=(100,'K'), Tmax=(1191.4,'K')), NASAPolynomial(coeffs=[14.0127,0.0122839,-5.66398e-06,1.10514e-09,-7.90106e-14,-2558.29,-67.4262], Tmin=(1191.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.16824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(12dioxetane) - ring(Cyclobutane) + ring(Ring) + ring(12dioxetane)"""),
)

species(
    label = '[O]C1C(=O)[CH]C1(F)F(1458)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {10,S}
8  C u0 p0 c0 {4,D} {5,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-324.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76843,0.0330212,2.0685e-05,-5.67295e-08,2.51147e-11,-38896.8,24.4092], Tmin=(100,'K'), Tmax=(977.3,'K')), NASAPolynomial(coeffs=[16.982,0.00982564,-3.68272e-06,8.01114e-10,-6.65757e-14,-43736.3,-58.1842], Tmin=(977.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-324.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsCsFF) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(CCJC=O)"""),
)

species(
    label = '[O]C1[C](F)C=C1F(1459)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5 C u1 p0 c0 {1,S} {4,S} {7,S}
6 C u0 p0 c0 {2,S} {4,S} {7,D}
7 C u0 p0 c0 {5,S} {6,D} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-76.9062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,180,864.045,864.046,864.046,864.048,864.055,864.056,864.058],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92628,0.0338811,1.66589e-06,-3.10557e-08,1.50484e-11,-9164.58,20.4671], Tmin=(100,'K'), Tmax=(991.386,'K')), NASAPolynomial(coeffs=[14.6396,0.0101379,-4.09666e-06,8.52018e-10,-6.69499e-14,-13039.3,-47.5857], Tmin=(991.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.9062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(CsCdCsF1s)"""),
)

species(
    label = 'O=C1C(F)=CC1(O)F(1460)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,D} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-576.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58202,0.059186,-7.7958e-05,5.99163e-08,-1.95523e-11,-69199,19.639], Tmin=(100,'K'), Tmax=(734.347,'K')), NASAPolynomial(coeffs=[7.20219,0.0285738,-1.54305e-05,3.15341e-09,-2.2865e-13,-70024.5,-5.73988], Tmin=(734.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-576.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(CdCCF) + longDistanceInteraction_cyclic(Cd(F)-CO) + longDistanceInteraction_cyclic(Cs(F)-CO) + ring(Cyclobutene)"""),
)

species(
    label = '[O]C1C2(F)[CH]C1(F)O2(1461)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-146.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54397,0.0391599,4.78248e-06,-3.76663e-08,1.69606e-11,-17519,15.6199], Tmin=(100,'K'), Tmax=(1037.66,'K')), NASAPolynomial(coeffs=[17.1603,0.0130963,-6.88403e-06,1.53041e-09,-1.20564e-13,-22597.6,-69.1521], Tmin=(1037.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + polycyclic(s3_4_4_ane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[O]C1(F)C2OC1[C]2F(1462)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-120.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87838,0.0452865,-3.32551e-05,1.10584e-08,-1.45198e-12,-14449.4,16.7692], Tmin=(100,'K'), Tmax=(1751.6,'K')), NASAPolynomial(coeffs=[14.393,0.0167077,-8.78119e-06,1.74348e-09,-1.22482e-13,-18833.5,-50.6211], Tmin=(1751.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-120.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + polycyclic(s3_4_4_ane) + radical(O2sj(Cs-F1sCsCs)) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C1(F)[CH]C2(F)OC12(1463)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-195.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35942,0.0490178,-3.37115e-05,6.76891e-09,7.20676e-13,-23419.6,19.2084], Tmin=(100,'K'), Tmax=(1205.61,'K')), NASAPolynomial(coeffs=[15.7499,0.0146761,-7.66052e-06,1.58508e-09,-1.16608e-13,-27863.6,-56.9482], Tmin=(1205.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFO) + group(Cs-CsCsHH) + polycyclic(s2_3_4_ane) + radical(O2sj(Cs-F1sCsCs)) + radical(CCJCO)"""),
)

species(
    label = '[O]C1[C](F)C2OC12F(1464)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-177.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64574,0.0514403,-4.2974e-05,1.67409e-08,-2.60613e-12,-21207.1,17.8473], Tmin=(100,'K'), Tmax=(1494.28,'K')), NASAPolynomial(coeffs=[13.906,0.0186211,-1.00291e-05,2.04258e-09,-1.47029e-13,-24871.1,-46.2254], Tmin=(1494.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + polycyclic(s2_3_4_ane) + radical(CC(C)OJ) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C=C(F)[CH]C(=O)F(1466)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u1 p0 c0 {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {10,S}
8  C u0 p0 c0 {2,S} {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-434.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,271,519,563,612,1379,3010,987.5,1337.5,450,1655,611,648,830,1210,1753,449.362,449.442,449.726],'cm^-1')),
        HinderedRotor(inertia=(0.224995,'amu*angstrom^2'), symmetry=1, barrier=(32.2118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224228,'amu*angstrom^2'), symmetry=1, barrier=(32.2083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20452,0.0491706,-2.37044e-05,-9.02686e-09,7.08405e-12,-52091.6,24.9913], Tmin=(100,'K'), Tmax=(1088.74,'K')), NASAPolynomial(coeffs=[17.9219,0.0124025,-7.0108e-06,1.54803e-09,-1.1967e-13,-57192.7,-63.7916], Tmin=(1088.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-434.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(COCsFO) + radical(C=COJ) + radical(CCJC=O)"""),
)

species(
    label = '[O]C(F)(C=O)C=[C]F(1467)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {4,D} {5,S} {10,S}
8  C u1 p0 c0 {2,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-201.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.960446,'amu*angstrom^2'), symmetry=1, barrier=(22.0826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966163,'amu*angstrom^2'), symmetry=1, barrier=(22.214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881316,0.0757134,-0.000126388,1.11102e-07,-3.82099e-11,-24122.7,25.6448], Tmin=(100,'K'), Tmax=(820.557,'K')), NASAPolynomial(coeffs=[9.16995,0.0243935,-1.26205e-05,2.46025e-09,-1.70848e-13,-25115.5,-10.4647], Tmin=(820.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C1C(=O)C=C1F(1468)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u0 p0 c0 {5,S} {6,D} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-102.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,323,467,575,827,1418,180,180,244.046,1073.13,1073.16,1073.17,1073.17,1073.22,1073.24,1073.27,2166.93,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20162,0.0330707,-1.50227e-05,-2.52963e-09,2.48237e-12,-12308.7,20.5136], Tmin=(100,'K'), Tmax=(1188.01,'K')), NASAPolynomial(coeffs=[10.9638,0.0152845,-7.35799e-06,1.47016e-09,-1.06157e-13,-15217.3,-26.7477], Tmin=(1188.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCsCdF) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C1(F)C=C(F)C1=O(1469)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {8,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {5,S} {8,D} {9,S}
7 C u0 p0 c0 {4,D} {5,S} {8,S}
8 C u0 p0 c0 {2,S} {6,D} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-332.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,288,410,724,839,1320,180,180,180,180,1325.73,1325.74,1326.12,1326.4],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58326,0.0594176,-9.45266e-05,8.74893e-08,-3.29111e-11,-39885.3,20.3944], Tmin=(100,'K'), Tmax=(753.687,'K')), NASAPolynomial(coeffs=[6.09278,0.0274704,-1.4995e-05,3.03216e-09,-2.16709e-13,-40337.5,1.4238], Tmin=(753.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCCF) + ring(Cyclobutene) + radical(C=OCOJ) + longDistanceInteraction_cyclic(Cd(F)-CO) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = 'O=C1[CH][C](F)C1=O(1470)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {7,D}
3 O u0 p2 c0 {6,D}
4 C u1 p0 c0 {5,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {4,S} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {7,S}
7 C u0 p0 c0 {2,D} {4,S} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-35.0852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,234,347,1316,1464,180,409.618,979.95,980.003,980.016,980.019,980.028,980.042,980.044,980.045,980.091,980.094],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48267,0.0248229,5.81468e-06,-2.38409e-08,9.8949e-12,-4157.83,20.4323], Tmin=(100,'K'), Tmax=(1064.61,'K')), NASAPolynomial(coeffs=[10.562,0.0141213,-6.79953e-06,1.39943e-09,-1.04399e-13,-6991.92,-24.2828], Tmin=(1064.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.0852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(CsCCFH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclobutane) + radical(CCJCC=O) + radical(CsCOCsF1s)"""),
)

species(
    label = '[O]C1(F)C=[C]C1=O(1471)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 C u0 p0 c0 {4,S} {7,D} {8,S}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (121.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,180,180,180,228.285,839.319,840.308,840.455,842.067,842.279,2683.91],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88138,0.05187,-8.27398e-05,7.57742e-08,-2.84367e-11,14741.1,19.0681], Tmin=(100,'K'), Tmax=(725.577,'K')), NASAPolynomial(coeffs=[6.34399,0.0226278,-1.26935e-05,2.60061e-09,-1.87448e-13,14215.7,-0.188132], Tmin=(725.577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(C=OCOJ) + radical(C=CJC=O) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = '[O]C1C(=O)[C]=C1F(1472)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (141.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,246,474,533,1155,180,180,858.59,858.593,858.598,858.685,858.774,858.84,858.897,859.239,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17838,0.0366011,-3.00203e-05,1.12316e-08,-1.64363e-12,17031.2,21.3012], Tmin=(100,'K'), Tmax=(1617.68,'K')), NASAPolynomial(coeffs=[12.8894,0.0101162,-5.46197e-06,1.11082e-09,-7.95353e-14,13565.8,-35.5253], Tmin=(1617.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCsCdF) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(C=OCOJ) + radical(C=CJC=O)"""),
)

species(
    label = '[O]C1C#CC1([O])F(1473)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 C u0 p0 c0 {4,S} {7,T}
7 C u0 p0 c0 {5,S} {6,T}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (47.5019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,360.685,360.685,360.686,360.686,360.687,360.687,360.687,360.687,360.687,360.687],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940216,0.0566002,-6.33173e-05,3.32507e-08,-6.53377e-12,5832.48,3.80073], Tmin=(100,'K'), Tmax=(1394.32,'K')), NASAPolynomial(coeffs=[17.0497,0.00421183,-3.16628e-07,-4.73852e-11,5.90269e-15,1940.25,-77.1211], Tmin=(1394.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.5019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + ring(Ring) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1(F)[CH]C(F)=C1O(1474)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {8,D}
7  C u1 p0 c0 {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {6,D} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-319.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431601,0.0686822,-7.28955e-05,3.60882e-08,-6.83038e-12,-38237.8,20.2842], Tmin=(100,'K'), Tmax=(1304.63,'K')), NASAPolynomial(coeffs=[20.0698,0.00847095,-3.66703e-06,7.12128e-10,-5.13774e-14,-43361.8,-79.6808], Tmin=(1304.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(C=CCJCO) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C1=C(F)[CH]C1(O)F(1475)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-411.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04476,0.0639499,-7.01008e-05,3.78e-08,-8.03355e-12,-49399.1,19.7722], Tmin=(100,'K'), Tmax=(1143.36,'K')), NASAPolynomial(coeffs=[14.5789,0.0166002,-7.9804e-06,1.57835e-09,-1.13407e-13,-52493.9,-47.3353], Tmin=(1143.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-411.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(C=C(C)OJ) + radical(C=CCJCO) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C1C(F)=[C]C1(O)F(1476)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {4,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-179.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2950,1000,246,474,533,1155,468.165,468.413,468.414,468.458,468.571,468.622,468.71,468.722],'cm^-1')),
        HinderedRotor(inertia=(0.000766846,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02429,0.0601432,-6.20869e-05,3.06401e-08,-5.87347e-12,-21500.1,23.7106], Tmin=(100,'K'), Tmax=(1274.19,'K')), NASAPolynomial(coeffs=[16.2966,0.0122,-5.64796e-06,1.11115e-09,-7.98628e-14,-25392.1,-53.6704], Tmin=(1274.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[O]C1(F)[C]=C(F)C1O(1477)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-179.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,1000,1380,1390,370,380,2900,435,246,474,533,1155,468.165,468.413,468.415,468.458,468.571,468.622,468.711,468.722],'cm^-1')),
        HinderedRotor(inertia=(0.00076685,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02429,0.0601432,-6.20869e-05,3.06401e-08,-5.87347e-12,-21500.1,23.7106], Tmin=(100,'K'), Tmax=(1274.19,'K')), NASAPolynomial(coeffs=[16.2966,0.0122,-5.64796e-06,1.11115e-09,-7.98628e-14,-25392.1,-53.6704], Tmin=(1274.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[O]C1[C](F)C=C1OF(1478)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {6,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {5,S} {8,D}
7  C u1 p0 c0 {1,S} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (74.8399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,3150,900,1100,346,659,817,1284,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24795,0.0506884,-3.32672e-05,1.36222e-09,4.04624e-12,9108.78,25.3301], Tmin=(100,'K'), Tmax=(1044.17,'K')), NASAPolynomial(coeffs=[16.2631,0.012249,-5.45795e-06,1.10793e-09,-8.30747e-14,4932.92,-52.7394], Tmin=(1044.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(CsCdCsF1s)"""),
)

species(
    label = '[O]C1=C[C](F)C1OF(1479)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {1,S} {5,S} {8,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-41.0145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,346,659,817,1284,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14644,0.0508334,-2.64404e-05,-1.22433e-08,1.07001e-11,-4819.44,25.0892], Tmin=(100,'K'), Tmax=(974.115,'K')), NASAPolynomial(coeffs=[18.3585,0.00775553,-2.60664e-06,5.3156e-10,-4.32657e-14,-9482.22,-64.2191], Tmin=(974.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.0145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(C=C(C)OJ) + radical(CsCdCsF1s)"""),
)

species(
    label = '[O]C1[C]=CC1(F)OF(1480)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (163.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93121,0.059544,-5.78592e-05,2.64721e-08,-4.68712e-12,19759.5,24.1242], Tmin=(100,'K'), Tmax=(1378.26,'K')), NASAPolynomial(coeffs=[17.5749,0.0112411,-5.2906e-06,1.04496e-09,-7.5004e-14,15171.6,-61.512], Tmin=(1378.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[O]C1(F)C=[C]C1OF(1481)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (163.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,3150,900,1100,1380,1390,370,380,2900,435,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931211,0.059544,-5.78592e-05,2.64721e-08,-4.68712e-12,19759.5,24.1242], Tmin=(100,'K'), Tmax=(1378.27,'K')), NASAPolynomial(coeffs=[17.5749,0.0112411,-5.29059e-06,1.04496e-09,-7.50038e-14,15171.6,-61.5121], Tmin=(1378.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(cyclobutene-vinyl)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C(F)C([O])F(378)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5 C u0 p0 c0 {2,S} {4,S} {6,D}
6 C u1 p0 c0 {5,D} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-79.3897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,246,474,533,1155,3120,650,792.5,1650,180,1192.01],'cm^-1')),
        HinderedRotor(inertia=(0.320481,'amu*angstrom^2'), symmetry=1, barrier=(7.36848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3380.05,'J/mol'), sigma=(5.44545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.96 K, Pc=47.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35812,0.0392947,-4.19267e-05,2.43185e-08,-5.88015e-12,-9491.9,18.0254], Tmin=(100,'K'), Tmax=(980.032,'K')), NASAPolynomial(coeffs=[7.66394,0.0176388,-8.78068e-06,1.77069e-09,-1.28302e-13,-10531.9,-7.46499], Tmin=(980.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.3897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCdH)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=C(F)O[C](F)C=O(1482)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u1 p0 c0 {1,S} {3,S} {7,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u0 p0 c0 {4,D} {5,S} {9,S}
8  C u1 p0 c0 {6,D} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-264.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17744,0.0644597,-7.52223e-05,4.44112e-08,-1.04869e-11,-31685.1,22.8959], Tmin=(100,'K'), Tmax=(1025.07,'K')), NASAPolynomial(coeffs=[12.5962,0.0199013,-1.00187e-05,2.00479e-09,-1.44487e-13,-34026.1,-32.4756], Tmin=(1025.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCFO) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCOF1sO2s) + radical(Cdj(Cd-F1sO2s)(H))"""),
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
    label = 'O=CC(=O)F(335)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C(F)=CC(=O)F(1483)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,D} {6,S} {8,S}
5 C u0 p0 c0 {1,S} {4,D} {7,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 C u2 p0 c0 {5,S} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-193.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,323,467,575,827,1418,255,533,799,832,1228,426.549,427.177,429.226,430.684],'cm^-1')),
        HinderedRotor(inertia=(0.413992,'amu*angstrom^2'), symmetry=1, barrier=(52.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400648,'amu*angstrom^2'), symmetry=1, barrier=(52.591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74072,0.0537618,-6.09553e-05,4.13702e-08,-1.20526e-11,-23242.3,20.9481], Tmin=(100,'K'), Tmax=(816.112,'K')), NASAPolynomial(coeffs=[7.00485,0.0279607,-1.35332e-05,2.63175e-09,-1.858e-13,-24101.5,-3.37853], Tmin=(816.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + group(Cd-Cd(CO)H) + group(COCFO) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=C(F)C1OC=C1F(1484)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {10,S}
8  C u0 p0 c0 {2,S} {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-510.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65808,0.0385701,-8.23398e-07,-3.19113e-08,1.58031e-11,-61328.7,22.526], Tmin=(100,'K'), Tmax=(995.938,'K')), NASAPolynomial(coeffs=[16.0599,0.0108739,-4.51293e-06,9.50773e-10,-7.50305e-14,-65692.4,-54.4016], Tmin=(995.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(CdCsCdF) + group(Cds-CdsOsH) + group(COCsFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cd-Cd(F)-Cs-O2s)"""),
)

species(
    label = 'C=C(F)C(=O)C(=O)F(1485)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,D}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,D} {6,S} {7,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {4,D} {5,S}
8  C u0 p0 c0 {6,D} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-643.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71756,0.050645,-4.74555e-05,2.284e-08,-4.48741e-12,-77304.4,22.2006], Tmin=(100,'K'), Tmax=(1203.66,'K')), NASAPolynomial(coeffs=[10.8565,0.0202745,-9.60768e-06,1.87728e-09,-1.3345e-13,-79504.4,-23.5835], Tmin=(1203.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(COCFO) + group(Cds-CdsHH)"""),
)

species(
    label = '[O]C1[C](F)OC=C1F(1486)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-272.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39936,0.0457869,-1.55426e-05,-2.41514e-08,1.60048e-11,-32691.2,19.3594], Tmin=(100,'K'), Tmax=(922.62,'K')), NASAPolynomial(coeffs=[17.1918,0.00700959,-7.69327e-07,5.34827e-11,-5.28674e-15,-36868.9,-62.4065], Tmin=(922.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[CH]=C(F)C1OC1([O])F(1487)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-107.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38333,0.0528058,-4.84284e-05,2.16186e-08,-3.80844e-12,-12868.1,24.8503], Tmin=(100,'K'), Tmax=(1365.64,'K')), NASAPolynomial(coeffs=[14.3028,0.0149639,-6.86326e-06,1.32755e-09,-9.38598e-14,-16396.7,-41.5048], Tmin=(1365.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFOO) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=C(F)C(=O)[C](O)F(1488)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {5,D}
5  C u0 p0 c0 {4,D} {6,S} {7,S}
6  C u1 p0 c0 {1,S} {3,S} {5,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {10,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-300.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,280,501,1494,1531,244,384,691,1241,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00258098,'amu*angstrom^2'), symmetry=1, barrier=(3.75522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25214,'amu*angstrom^2'), symmetry=1, barrier=(51.7811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25443,'amu*angstrom^2'), symmetry=1, barrier=(51.8337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09066,0.0714024,-0.000117202,1.0606e-07,-3.79605e-11,-35997.4,24.3208], Tmin=(100,'K'), Tmax=(803.901,'K')), NASAPolynomial(coeffs=[7.36369,0.0283779,-1.48832e-05,2.93129e-09,-2.05367e-13,-36624.3,-2.19986], Tmin=(803.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-Cds)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-CdsHH) + radical(CsCOF1sO2s) + radical(Cdj(Cd-F1sCO)(H))"""),
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
    label = 'CFO(51)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=C(F)C=O(471)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u0 p0 c0 {2,D} {3,S} {6,S}
5 C u1 p0 c0 {3,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-6.6029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([244,384,691,1241,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.03426,'amu*angstrom^2'), symmetry=1, barrier=(23.7797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0457,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88566,0.00800032,9.79848e-05,-2.62824e-07,2.05152e-10,-790.516,10.04], Tmin=(10,'K'), Tmax=(444.064,'K')), NASAPolynomial(coeffs=[4.44622,0.0222292,-1.51983e-05,4.85875e-09,-5.87968e-13,-1030.38,5.65049], Tmin=(444.064,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-6.6029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[CH]DC(F)CDO""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=C(F)C(=O)C(=O)F(1489)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {3,D} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {5,S} {8,D}
7 C u0 p0 c0 {2,S} {4,D} {5,S}
8 C u1 p0 c0 {6,D} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-384.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,244,384,691,1241,286,619,818,1246,1924,3120,650,792.5,1650,180,1682.56],'cm^-1')),
        HinderedRotor(inertia=(0.945892,'amu*angstrom^2'), symmetry=1, barrier=(21.7479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.48771,'amu*angstrom^2'), symmetry=1, barrier=(57.1974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3752,0.0623196,-8.69006e-05,6.37286e-08,-1.87905e-11,-46213.9,22.9791], Tmin=(100,'K'), Tmax=(826.633,'K')), NASAPolynomial(coeffs=[10.0791,0.0202022,-1.04743e-05,2.09154e-09,-1.49405e-13,-47652.9,-17.3548], Tmin=(826.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-384.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(COCFO) + group(Cds-CdsHH) + radical(Cdj(Cd-F1sCO)(H))"""),
)

species(
    label = 'C#CC([O])C(=O)F(1490)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 C u0 p0 c0 {4,S} {7,T}
7 C u0 p0 c0 {6,T} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-114.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,2175,525,750,770,3400,2100,439.969,439.979],'cm^-1')),
        HinderedRotor(inertia=(0.0592916,'amu*angstrom^2'), symmetry=1, barrier=(8.14859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63354,'amu*angstrom^2'), symmetry=1, barrier=(87.0402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66497,0.0493925,-5.6208e-05,3.28494e-08,-7.52351e-12,-13703.1,23.9586], Tmin=(100,'K'), Tmax=(1072.44,'K')), NASAPolynomial(coeffs=[11.5462,0.0125373,-4.65927e-06,8.04855e-10,-5.35169e-14,-15822.5,-24.4039], Tmin=(1072.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(COCsFO) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(509)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.101,381.695],'cm^-1')),
        HinderedRotor(inertia=(0.482775,'amu*angstrom^2'), symmetry=1, barrier=(49.7784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80088,0.0290439,-3.02832e-05,1.66615e-08,-3.81631e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.22,'K')), NASAPolynomial(coeffs=[6.95396,0.0128879,-6.71487e-06,1.38086e-09,-1.01081e-13,-26608.7,-6.32755], Tmin=(1028.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
)

species(
    label = 'C#CC([O])=C([O])F(1491)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {2,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {4,D}
6 C u0 p0 c0 {4,S} {7,T}
7 C u0 p0 c0 {6,T} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-31.1229,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,326,540,652,719,1357,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.85925,'amu*angstrom^2'), symmetry=1, barrier=(42.7477,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60084,0.0567468,-8.0948e-05,5.9597e-08,-1.74934e-11,-3660.37,19.5408], Tmin=(100,'K'), Tmax=(832.896,'K')), NASAPolynomial(coeffs=[9.95947,0.0166036,-8.65129e-06,1.7283e-09,-1.23387e-13,-5052.73,-19.2563], Tmin=(832.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.1229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C(F)C(=O)[C]=O(1492)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {4,S} {7,D}
6 C u1 p0 c0 {3,D} {4,S}
7 C u1 p0 c0 {5,D} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (38.4605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,244,384,691,1241,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0388648,'amu*angstrom^2'), symmetry=1, barrier=(12.2424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.47504,'amu*angstrom^2'), symmetry=1, barrier=(56.9061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65424,0.0561439,-8.28377e-05,6.18748e-08,-1.74729e-11,4705.99,21.791], Tmin=(100,'K'), Tmax=(673.263,'K')), NASAPolynomial(coeffs=[9.33947,0.0166017,-8.36869e-06,1.6312e-09,-1.14176e-13,3532.5,-13.275], Tmin=(673.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.4605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(Cdj(Cd-F1sCO)(H)) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C(F)=C(O)C(=O)F(1493)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {9,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {1,S} {5,D} {8,S}
7  C u0 p0 c0 {2,S} {4,D} {5,S}
8  C u2 p0 c0 {6,S} {10,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-373.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732245,0.0751521,-9.58188e-05,6.36342e-08,-1.68924e-11,-44825.7,23.8397], Tmin=(100,'K'), Tmax=(918.106,'K')), NASAPolynomial(coeffs=[12.6308,0.0233165,-1.1136e-05,2.14784e-09,-1.50934e-13,-47010.7,-32.5481], Tmin=(918.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-373.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COCFO) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(F)C([O])=C([O])F(1494)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {4,S} {5,D}
8  C u0 p0 c0 {6,D} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11175,0.0690044,-9.71874e-05,7.31911e-08,-2.22966e-11,-48370.6,22.9756], Tmin=(100,'K'), Tmax=(799.074,'K')), NASAPolynomial(coeffs=[10.1304,0.0238601,-1.24457e-05,2.49283e-09,-1.78293e-13,-49812,-18.5114], Tmin=(799.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C(OF)C(=O)F(1495)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u1 p0 c0 {5,S} {8,D}
8  C u1 p0 c0 {7,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (95.1154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,1685,370,3120,650,792.5,1650,213.665,805.714],'cm^-1')),
        HinderedRotor(inertia=(0.586953,'amu*angstrom^2'), symmetry=1, barrier=(18.7452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583162,'amu*angstrom^2'), symmetry=1, barrier=(18.6964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224507,'amu*angstrom^2'), symmetry=1, barrier=(7.17694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20949,0.0656533,-8.56399e-05,5.72896e-08,-1.53445e-11,11536.6,28.2593], Tmin=(100,'K'), Tmax=(907.871,'K')), NASAPolynomial(coeffs=[11.5299,0.0201826,-1.05127e-05,2.12249e-09,-1.53202e-13,9662.67,-20.5332], Tmin=(907.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.1154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[O]C([C]=CF)C(=O)F(1496)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u0 p0 c0 {2,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-234.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49889,0.05923,-7.72596e-05,5.47755e-08,-1.58093e-11,-28123.1,28.2787], Tmin=(100,'K'), Tmax=(840.197,'K')), NASAPolynomial(coeffs=[9.30143,0.0220845,-1.09453e-05,2.15857e-09,-1.53512e-13,-29434.2,-8.00564], Tmin=(840.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-234.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(F)C([C]=O)OF(1497)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u1 p0 c0 {4,D} {5,S}
8  C u1 p0 c0 {6,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (73.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,246,474,533,1155,1855,455,950,3120,650,792.5,1650,292.327,292.371],'cm^-1')),
        HinderedRotor(inertia=(0.392771,'amu*angstrom^2'), symmetry=1, barrier=(23.8422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210035,'amu*angstrom^2'), symmetry=1, barrier=(23.8474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233214,'amu*angstrom^2'), symmetry=1, barrier=(14.1655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02745,0.0676511,-8.77572e-05,5.64442e-08,-1.4262e-11,8992.47,28.6464], Tmin=(100,'K'), Tmax=(969.263,'K')), NASAPolynomial(coeffs=[13.5149,0.0161178,-8.00686e-06,1.59197e-09,-1.14254e-13,6571.71,-31.2087], Tmin=(969.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CdCsCdF) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[O]C([C]=O)C(F)=CF(1498)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {10,S}
8  C u1 p0 c0 {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-242.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24248,0.0642351,-8.83157e-05,6.3132e-08,-1.79397e-11,-29011.3,28.6304], Tmin=(100,'K'), Tmax=(861.282,'K')), NASAPolynomial(coeffs=[11.0471,0.0187007,-9.0149e-06,1.75108e-09,-1.23281e-13,-30700.3,-17.2073], Tmin=(861.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    E0 = (93.4005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (304.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (273.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (164.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (150.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-90.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (45.1191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (299.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (114.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (56.2649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (74.2211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (232.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (408.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (439.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (128.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (383.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (508.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (177.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (46.4436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (362.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (161.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (359.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (185.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (127.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (136.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-91.4706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-25.3289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-143.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (113.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (131.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (234.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-81.9747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (93.5164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (123.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (305.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (335.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (199.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (265.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (159.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (235.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (311.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (249.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (285.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (152.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (98.0161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (216.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (150.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-88.2707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-20.8135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-5.33793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (9.93146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (-1.36656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (-127.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (-151.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (-97.7446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (35.8382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (-54.8775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (58.0447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (133.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (144.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (196.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (-26.2108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (-40.6264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (12.4311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (-28.3684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (213.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (164.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (276.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (299.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (220.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (-35.4985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (-180.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (99.2842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (-172.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (-117.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (123.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (-132.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (-57.6624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (7.01457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (-27.8321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (-66.4991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (-100.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (70.4574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (-60.4105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (186.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (43.3904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (76.7853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (-55.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (11.4378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (242.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (-0.414337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (227.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (42.0684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', 'FC1=CC(F)=C1(307)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.23921e+14,'m^3/(mol*s)'), n=-2.03695, Ea=(22.8417,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.8381518263401699, var=29.6634264973539, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_N-4R!H-u0_Ext-2CNS-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_N-4R!H-u0_Ext-2CNS-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]O[C](F)C1C=C1F(936)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OC1=CC(F)(F)[CH]1(937)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(253.187,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['O=O(207)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(264.03,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 264.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', '[O]C1(F)[CH]C(F)=C1(938)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['O(6)', 'FC1=CC2(F)OC12(940)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.54e+09,'s^-1'), n=1.06, Ea=(144.278,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOJ] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['HO2(11)', 'FC1=[C]C(F)=C1(576)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.26e+09,'s^-1'), n=1.11, Ea=(398.636,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['[O]OC1(F)C2[C](F)C21(941)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['FC12[CH]C(F)([CH]1)OO2(942)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.17321e+11,'s^-1'), n=0.14784, Ea=(155.424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['F[C]1[CH]C2(F)OOC12(943)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.21728e+11,'s^-1'), n=0.310877, Ea=(173.38,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 172.0 to 173.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(F)C=C(F)O[O](944)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]OC1=CC(F)=C1(945)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(7.52229,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]OC1(F)C=C=C1(946)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(5.27814,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.4916e+06,'m^3/(mol*s)'), n=-0.00476945, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_2C-inRing_Ext-2C-R_Ext-4R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_2C-inRing_Ext-2C-R_Ext-4R!H-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', '[O]OC1=C=C(F)[CH]1(947)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(119.434,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', '[O]OC1(F)[C]=C=C1(948)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(190.596,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC1(F)[C]=C(F)C1(949)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['OOC1(F)[C]=C(F)[CH]1(950)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(548,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['FOO[C]1[CH]C(F)=C1(951)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(82.4259,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]O[C]1C=C(F)C1F(952)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(130.828,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['FOOC1(F)[CH][C]=C1(953)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(54.2652,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]OC1(F)C=[C]C1F(954)'],
    products = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(114.826,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FC1=CC2(OO2)C1F(961)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC1(F)C=C2OOC21(962)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(283.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]1C=C(F)[CH]OO1(963)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC1[C](F)C=C1F(580)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[C]=CC1(F)[CH]OO1(965)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;Cdsinglepri_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(F)C1OO[C]1F(966)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_noH;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=O(207)', 'FC1=CC(F)=C1(307)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.83584e-11,'m^3/(mol*s)'), n=4.15719, Ea=(69.212,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00023424662986478353, var=0.0043250921906934775, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-5COCSCdCddCtN3dN3tN5dcN5tcO2dS2dS4dS4tS6dS6tS6tdS6tt->Cd',), comment="""Estimated from node Root_N-5COCSCdCddCtN3dN3tN5dcN5tcO2dS2dS4dS4tS6dS6tS6tdS6tt->Cd
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['FC1=COOC(F)=C1(967)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FC1[CH]C2(F)OO[C]12(968)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['F[C]1CC2(F)OO[C]12(969)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F(37)', 'FC1=C[C]2OOC21(970)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction36',
    reactants = ['F(37)', 'FC12C=[C]C1OO2(971)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(5)', 'FC1=CC2(F)OO[C]12(972)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(5)', 'FC1=[C]C2(F)OOC12(973)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction39',
    reactants = ['FC1[C]C2(F)OOC12(974)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.03468e+20,'s^-1'), n=-2.18552, Ea=(34.6055,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5822770821690705, var=6.058307656543806, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C"""),
)

reaction(
    label = 'reaction40',
    reactants = ['FC1[C]C2OOC12F(975)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(84.6696,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction41',
    reactants = ['HF(38)', 'FC1=CC2=C1OO2(976)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(106.26,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction42',
    reactants = ['HF(38)', 'FC12C=C=C1OO2(977)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(134.648,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction43',
    reactants = ['HF(38)', 'FC1=C=C2OOC12(978)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(119.023,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction44',
    reactants = ['HF(38)', 'FC12C#CC1OO2(979)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(376.379,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1C(=O)[CH]C1(F)F(1458)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(249.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction46',
    reactants = ['O(6)', '[O]C1[C](F)C=C1F(1459)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction47',
    reactants = ['O(6)', '[O]C1(F)[CH]C(F)=C1(938)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['O=C1C(F)=CC1(O)F(1460)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1C2(F)[CH]C1(F)O2(1461)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c4_beta;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1(F)C2OC1[C]2F(1462)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.44078e+13,'s^-1'), n=-0.00629116, Ea=(146.333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_cyclic;multiplebond_intra;radadd_intra_O] + [R4_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c4_beta;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1(F)[CH]C2(F)OC12(1463)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.22471e+12,'s^-1'), n=0.01909, Ea=(161.602,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c4_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1[C](F)C2OC12F(1464)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.91816e+12,'s^-1'), n=-0.119691, Ea=(150.304,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_DS_CO;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[O]C=C(F)[CH]C(=O)F(1466)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(232.64,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[O]C(F)(C=O)C=[C]F(1467)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_H;radadd_intra_cdsingle]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['F(37)', '[O]C1C(=O)C=C1F(1468)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(15.7254,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction57',
    reactants = ['H(5)', '[O]C1(F)C=C(F)C1=O(1469)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(15.4723,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction58',
    reactants = ['HF(38)', 'O=C1[CH][C](F)C1=O(1470)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(324.099,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction59',
    reactants = ['HF(38)', '[O]C1(F)C=[C]C1=O(1471)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(242.654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction60',
    reactants = ['HF(38)', '[O]C1C(=O)[C]=C1F(1472)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(234.682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction61',
    reactants = ['HF(38)', '[O]C1C#CC1([O])F(1473)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(379.495,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1(F)[CH]C(F)=C1O(1474)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1=C(F)[CH]C1(O)F(1475)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.99911e+07,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;O_rad_out;XH_out] for rate rule [R3H_SS_23cy4;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[O]C1C(F)=[C]C1(O)F(1476)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_12cy4;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[O]C1(F)[C]=C(F)C1O(1477)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[O]C1[C](F)C=C1OF(1478)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(88.3018,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[O]C1=C[C](F)C1OF(1479)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(155.404,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[O]C1[C]=CC1(F)OF(1480)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(63.3302,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[O]C1(F)C=[C]C1OF(1481)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(85.9073,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction70',
    reactants = ['CO(13)', '[CH]=C(F)C([O])F(378)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(368.875,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[CH]=C(F)O[C](F)C=O(1482)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(145.193,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['C2HF(58)', 'O=CC(=O)F(335)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction73',
    reactants = ['O(6)', '[CH]C(F)=CC(=O)F(1483)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['O=C(F)C1OC=C1F(1484)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['C=C(F)C(=O)C(=O)F(1485)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[CH]=C(F)C1OO[C]1F(966)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 303.2 to 303.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction77',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[O]C1[C](F)OC=C1F(1486)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(3.00067e+09,'s^-1'), n=0.569069, Ea=(47.8102,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_DS_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction78',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[CH]=C(F)C1OC1([O])F(1487)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(123.029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 123.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH]=C(F)C(=O)[C](O)F(1488)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(256.981,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction80',
    reactants = ['[CH]=[C]F(252)', 'O=CC(=O)F(335)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(54.6935,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction81',
    reactants = ['CFO(51)', '[CH]=C(F)C=O(471)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(80.3185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction82',
    reactants = ['H(5)', '[CH]=C(F)C(=O)C(=O)F(1489)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(22.4406,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction83',
    reactants = ['F(37)', 'C#CC([O])C(=O)F(1490)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.0689,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction84',
    reactants = ['C2HF(58)', '[O][CH]C(=O)F(509)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(8.59146,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction85',
    reactants = ['[CH]=[C]F(252)', '[O][CH]C(=O)F(509)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction86',
    reactants = ['HF(38)', 'C#CC([O])=C([O])F(1491)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(305.483,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction87',
    reactants = ['HF(38)', '[CH]=C(F)C(=O)[C]=O(1492)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(269.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction88',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[CH]C(F)=C(O)C(=O)F(1493)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction89',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['C=C(F)C([O])=C([O])F(1494)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction90',
    reactants = ['[CH]=[C]C(OF)C(=O)F(1495)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(97.2399,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction91',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[O]C([C]=CF)C(=O)F(1496)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(180.277,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction92',
    reactants = ['[CH]=C(F)C([C]=O)OF(1497)'],
    products = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction93',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[O]C([C]=O)C(F)=CF(1498)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(222.759,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

network(
    label = 'PDepNetwork #121',
    isomers = [
        '[O]OC1(F)[CH]C(F)=C1(581)',
        'FC1=CC2(F)OOC12(939)',
        '[O]C1C(F)=CC1([O])F(964)',
        '[CH]=C(F)C([O])C(=O)F(1465)',
    ],
    reactants = [
        ('O2(2)', 'FC1=CC(F)=C1(307)'),
        ('O=O(207)', 'FC1=CC(F)=C1(307)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #121',
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

