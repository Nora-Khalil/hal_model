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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14284,0.052484,-3.48088e-05,1.04848e-09,4.54051e-12,-24160.8,22.6183], Tmin=(100,'K'), Tmax=(1035.18,'K')), NASAPolynomial(coeffs=[17.0604,0.0114786,-5.098e-06,1.04607e-09,-7.9314e-14,-28554.8,-60.0307], Tmin=(1035.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    E0 = (76.0776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (85.5742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-142.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-142.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-76.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-194.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (62.3751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (80.2309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (182.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-133.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (42.2653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (72.5787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (254.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (284.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (148.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (214.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (108.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (183.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (259.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (198.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (233.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (101.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction24',
    reactants = ['FC1=CC2(OO2)C1F(961)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction25',
    reactants = ['FC1(F)C=C2OOC21(962)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.7779e+11,'s^-1'), n=0.725184, Ea=(283.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]1C=C(F)[CH]OO1(963)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC1(F)[CH]C(F)=C1(581)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC1[C](F)C=C1F(580)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['F[C]=CC1(F)[CH]OO1(965)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;Cdsinglepri_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(F)C1OO[C]1F(966)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_noH;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=O(207)', 'FC1=CC(F)=C1(307)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.83584e-11,'m^3/(mol*s)'), n=4.15719, Ea=(69.212,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00023424662986478353, var=0.0043250921906934775, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-5COCSCdCddCtN3dN3tN5dcN5tcO2dS2dS4dS4tS6dS6tS6tdS6tt->Cd',), comment="""Estimated from node Root_N-5COCSCdCddCtN3dN3tN5dcN5tcO2dS2dS4dS4tS6dS6tS6tdS6tt->Cd
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction32',
    reactants = ['FC1=COOC(F)=C1(967)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1,3-butadiene_backbone;C=C_1;C=C_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction33',
    reactants = ['FC1[CH]C2(F)OO[C]12(968)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['F[C]1CC2(F)OO[C]12(969)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F(37)', 'FC1=C[C]2OOC21(970)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction36',
    reactants = ['F(37)', 'FC12C=[C]C1OO2(971)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(5)', 'FC1=CC2(F)OO[C]12(972)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(5)', 'FC1=[C]C2(F)OOC12(973)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction39',
    reactants = ['FC1[C]C2(F)OOC12(974)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.03468e+20,'s^-1'), n=-2.18552, Ea=(34.6055,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.5822770821690705, var=6.058307656543806, Tref=1000.0, N=2, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing_Ext-4CCl-R_Ext-5R!H-R_Sp-5R!H-1C"""),
)

reaction(
    label = 'reaction40',
    reactants = ['FC1[C]C2OOC12F(975)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(84.6696,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction41',
    reactants = ['HF(38)', 'FC1=CC2=C1OO2(976)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(106.26,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction42',
    reactants = ['HF(38)', 'FC12C=C=C1OO2(977)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(134.648,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction43',
    reactants = ['HF(38)', 'FC1=C=C2OOC12(978)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(119.023,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction44',
    reactants = ['HF(38)', 'FC12C#CC1OO2(979)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(376.379,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

network(
    label = 'PDepNetwork #185',
    isomers = [
        'FC1=CC2(F)OOC12(939)',
    ],
    reactants = [
        ('O=O(207)', 'FC1=CC(F)=C1(307)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #185',
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

