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
    E0 = (159.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (278.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (212.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-81.4479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-26.3321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (41.1251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (56.6007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (71.8701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (60.5721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-65.1972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-89.3148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-35.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (97.7769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (7.06114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (119.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (195.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (206.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (257.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (35.7278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (21.3123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (74.3698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (33.5703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (275.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (226.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (338.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (361.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction45',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1C(=O)[CH]C1(F)F(1458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(249.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction46',
    reactants = ['O(6)', '[O]C1[C](F)C=C1F(1459)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS2',
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
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['FC1=CC2(F)OOC12(939)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['O=C1C(F)=CC1(O)F(1460)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1C2(F)[CH]C1(F)O2(1461)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c4_beta;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1(F)C2OC1[C]2F(1462)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.44078e+13,'s^-1'), n=-0.00629116, Ea=(146.333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_cyclic;multiplebond_intra;radadd_intra_O] + [R4_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c4_beta;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1(F)[CH]C2(F)OC12(1463)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.22471e+12,'s^-1'), n=0.01909, Ea=(161.602,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c4_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1[C](F)C2OC12F(1464)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.91816e+12,'s^-1'), n=-0.119691, Ea=(150.304,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c4_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=C(F)C([O])C(=O)F(1465)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_DS_CO;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[O]C=C(F)[CH]C(=O)F(1466)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(232.64,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[O]C(F)(C=O)C=[C]F(1467)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_H;radadd_intra_cdsingle]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['F(37)', '[O]C1C(=O)C=C1F(1468)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(15.7254,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction57',
    reactants = ['H(5)', '[O]C1(F)C=C(F)C1=O(1469)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(15.4723,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction58',
    reactants = ['HF(38)', 'O=C1[CH][C](F)C1=O(1470)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(324.099,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction59',
    reactants = ['HF(38)', '[O]C1(F)C=[C]C1=O(1471)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(242.654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction60',
    reactants = ['HF(38)', '[O]C1C(=O)[C]=C1F(1472)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(234.682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction61',
    reactants = ['HF(38)', '[O]C1C#CC1([O])F(1473)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(379.495,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1(F)[CH]C(F)=C1O(1474)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[O]C1C(F)=CC1([O])F(964)'],
    products = ['[O]C1=C(F)[CH]C1(O)F(1475)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.99911e+07,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;O_rad_out;XH_out] for rate rule [R3H_SS_23cy4;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[O]C1C(F)=[C]C1(O)F(1476)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_12cy4;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[O]C1(F)[C]=C(F)C1O(1477)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[O]C1[C](F)C=C1OF(1478)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(88.3018,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[O]C1=C[C](F)C1OF(1479)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(155.404,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[O]C1[C]=CC1(F)OF(1480)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(63.3302,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[O]C1(F)C=[C]C1OF(1481)'],
    products = ['[O]C1C(F)=CC1([O])F(964)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(85.9073,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #251',
    isomers = [
        '[O]C1C(F)=CC1([O])F(964)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #251',
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

