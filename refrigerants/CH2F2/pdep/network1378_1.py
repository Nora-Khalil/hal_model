species(
    label = '[O]C(F)(F)O[CH]C(=O)F(3270)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u1 p0 c0 {4,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-769.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,3025,407.5,1350,352.5,611,648,830,1210,1753,180,180,180,741.86],'cm^-1')),
        HinderedRotor(inertia=(0.0887822,'amu*angstrom^2'), symmetry=1, barrier=(2.04128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794606,'amu*angstrom^2'), symmetry=1, barrier=(18.2695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.463936,'amu*angstrom^2'), symmetry=1, barrier=(18.2695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576174,0.0770895,-0.000102379,6.53001e-08,-1.6183e-11,-92409.9,28.0188], Tmin=(100,'K'), Tmax=(992.168,'K')), NASAPolynomial(coeffs=[16.1523,0.0142934,-7.44256e-06,1.50962e-09,-1.09643e-13,-95500.8,-47.0049], Tmin=(992.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-769.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(CsFFOO) + group(COCsFO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(CCsJOCs)"""),
)

species(
    label = 'CF2O(48)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-618.61,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CC(=O)F(711)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C([O])(F)F(150)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-372.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,464.782],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3646.55,'J/mol'), sigma=(6.01309,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.58 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92149,0.00463651,6.40183e-05,-1.39109e-07,8.57731e-11,-44748,10.2194], Tmin=(10,'K'), Tmax=(572.951,'K')), NASAPolynomial(coeffs=[4.84692,0.0160279,-1.25421e-05,4.35666e-09,-5.55356e-13,-45147.1,3.71309], Tmin=(572.951,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-372.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C([O])(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C(=O)F-2(3116)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-5.0725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([486,617,768,1157,1926,180,1655.08,1655.49],'cm^-1')),
        HinderedRotor(inertia=(0.0191603,'amu*angstrom^2'), symmetry=1, barrier=(5.31405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32766,0.0152939,-1.10759e-05,3.91583e-09,-5.61424e-13,-586.537,12.101], Tmin=(100,'K'), Tmax=(1580.39,'K')), NASAPolynomial(coeffs=[6.5355,0.00717482,-3.3698e-06,6.65151e-10,-4.72046e-14,-1600.47,-4.84309], Tmin=(1580.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
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
    label = 'O=C(F)[CH]O[C](F)F(2827)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {4,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 C u1 p0 c0 {2,S} {3,S} {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-610.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,493,600,700,1144,1293,243.578,244.083,246.058,2447.77],'cm^-1')),
        HinderedRotor(inertia=(0.234687,'amu*angstrom^2'), symmetry=1, barrier=(10.0373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241975,'amu*angstrom^2'), symmetry=1, barrier=(10.0163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644436,'amu*angstrom^2'), symmetry=1, barrier=(27.2967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888372,0.074163,-0.000124072,1.04806e-07,-3.46464e-11,-73370.2,26.6213], Tmin=(100,'K'), Tmax=(800.145,'K')), NASAPolynomial(coeffs=[11.1628,0.0179603,-9.6383e-06,1.90232e-09,-1.32983e-13,-74859.5,-19.6881], Tmin=(800.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFHO) + group(COCsFO) + radical(CCsJOCs) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O=C(F)C1OC(F)(F)O1(3354)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1111.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12621,0.0442466,-3.21239e-05,1.08472e-08,-1.4904e-12,-133630,22.8034], Tmin=(100,'K'), Tmax=(1616.49,'K')), NASAPolynomial(coeffs=[11.2267,0.0217274,-1.12274e-05,2.2291e-09,-1.57541e-13,-136573,-25.4714], Tmin=(1616.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1111.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(CsFFOO) + group(COCsFO) + ring(Cs-O2s-Cs(F)(F)-O2s)"""),
)

species(
    label = '[O]C(F)(F)OC1O[C]1F(3426)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
9  C u1 p0 c0 {3,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-619.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08752,0.0619772,-6.34043e-05,3.11918e-08,-6.04761e-12,-74374,25.6709], Tmin=(100,'K'), Tmax=(1245.24,'K')), NASAPolynomial(coeffs=[15.3365,0.0162061,-8.26918e-06,1.67409e-09,-1.21503e-13,-77922.7,-46.1976], Tmin=(1245.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFOO) + ring(Cs(O2)-O2s-Cs(F)) + radical(O2sj(Cs-F1sF1sO2s)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'F[C]1[CH]OC(F)(F)OO1(3281)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u1 p0 c0 {4,S} {9,S} {10,S}
9  C u1 p0 c0 {3,S} {6,S} {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-600.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.801123,0.0807981,-9.64801e-05,5.31689e-08,-1.04591e-11,-72004.7,26.9444], Tmin=(100,'K'), Tmax=(1536.92,'K')), NASAPolynomial(coeffs=[20.5504,0.000753013,5.5298e-06,-1.4415e-09,1.09497e-13,-75677.1,-75.8365], Tmin=(1536.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-600.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsFFOO) + ring(124trioxane) + radical(CCsJOCs) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)[CH]OC(F)(F)O1(3308)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
9  C u1 p0 c0 {5,S} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-788.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.123363,0.076181,-0.00010305,6.26526e-08,-1.35547e-11,-94732.5,22.6747], Tmin=(100,'K'), Tmax=(1384.5,'K')), NASAPolynomial(coeffs=[19.7205,-0.0031206,6.67067e-06,-1.64203e-09,1.24727e-13,-98121.6,-71.9118], Tmin=(1384.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-788.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsFFOO) + ring(1,3-Dioxolane) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOCs)"""),
)

species(
    label = '[O][C](F)F(237)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-255.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53788e-05,4.83907e-09,3.86776e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58999,0.00108969,-4.53602e-07,1.24872e-10,-1.13713e-14,-32130.4,-16.3888], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
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
    label = 'O=C(F)[CH]OC(=O)F(3427)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {3,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-756.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,482,664,788,1296,1923,314.813,314.815,314.818,314.819],'cm^-1')),
        HinderedRotor(inertia=(0.00170097,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334262,'amu*angstrom^2'), symmetry=1, barrier=(23.5091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334251,'amu*angstrom^2'), symmetry=1, barrier=(23.5091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27907,0.0603605,-7.18958e-05,4.10936e-08,-9.19401e-12,-90933.4,24.5936], Tmin=(100,'K'), Tmax=(1090.48,'K')), NASAPolynomial(coeffs=[14.1284,0.0132274,-7.06183e-06,1.45699e-09,-1.07003e-13,-93735.8,-38.5101], Tmin=(1090.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(COFOO) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][CH]C(=O)F(753)',
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
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,379.05,383.679],'cm^-1')),
        HinderedRotor(inertia=(0.48558,'amu*angstrom^2'), symmetry=1, barrier=(49.8051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8009,0.0290437,-3.02827e-05,1.6661e-08,-3.81614e-12,-25754.7,13.8242], Tmin=(100,'K'), Tmax=(1028.29,'K')), NASAPolynomial(coeffs=[6.95402,0.0128878,-6.7148e-06,1.38084e-09,-1.0108e-13,-26608.7,-6.32789], Tmin=(1028.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
)

species(
    label = '[O]C(F)(F)OC=C=O(3428)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
7 C u0 p0 c0 {3,S} {8,D} {9,S}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-568.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16553,'amu*angstrom^2'), symmetry=1, barrier=(26.7979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16493,'amu*angstrom^2'), symmetry=1, barrier=(26.784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0177191,0.0693183,-8.07821e-05,4.16678e-08,-7.8442e-12,-68237.1,25.8036], Tmin=(100,'K'), Tmax=(1510.07,'K')), NASAPolynomial(coeffs=[22.5735,-0.00264774,2.84096e-06,-6.0842e-10,4.19463e-14,-73656.2,-87.6988], Tmin=(1510.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-568.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsFFOO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(O2sj(Cs-F1sF1sO2s))"""),
)

species(
    label = 'O=C(F)[CH]O[C](F)OF(3429)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  C u1 p0 c0 {2,S} {4,S} {5,S}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-463.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,482,586,761,1411,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370867,0.0850924,-0.00012761,9.37154e-08,-2.69189e-11,-55565.9,29.5574], Tmin=(100,'K'), Tmax=(856.508,'K')), NASAPolynomial(coeffs=[14.783,0.0177873,-9.74183e-06,1.97397e-09,-1.41736e-13,-58034.8,-37.741], Tmin=(856.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-463.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(CsFHOO) + group(COCsFO) + radical(CCsJOCs) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O][C](F)OC(F)C(=O)F(3430)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {8,D}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-752.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,482,586,761,1411,180,180,328.661,1009.32,1009.97],'cm^-1')),
        HinderedRotor(inertia=(1.88431,'amu*angstrom^2'), symmetry=1, barrier=(43.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0175566,'amu*angstrom^2'), symmetry=1, barrier=(12.702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43454,'amu*angstrom^2'), symmetry=1, barrier=(32.983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2,0.0627678,-6.78277e-05,3.60943e-08,-7.66414e-12,-90403.4,28.3038], Tmin=(100,'K'), Tmax=(1133.45,'K')), NASAPolynomial(coeffs=[13.6078,0.0189797,-9.87834e-06,2.00968e-09,-1.46186e-13,-93216.1,-33.111], Tmin=(1133.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-752.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHOO) + group(COCsFO) + radical(O2sj(Cs-F1sO2sH)) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = '[O][C]=COC(F)(F)OF(3431)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {10,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-500.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,435,565,619,662,854,1178,1396,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24069,'amu*angstrom^2'), symmetry=1, barrier=(28.5259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24316,'amu*angstrom^2'), symmetry=1, barrier=(28.5828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23775,'amu*angstrom^2'), symmetry=1, barrier=(28.4583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.615942,0.0864865,-0.000106993,5.86874e-08,-1.18997e-11,-60071.5,31.0789], Tmin=(100,'K'), Tmax=(1339.51,'K')), NASAPolynomial(coeffs=[26.3666,-0.00262582,2.35729e-06,-4.93921e-10,3.3756e-14,-66534.2,-104.124], Tmin=(1339.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-500.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsFFOO) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)(F)OC(F)[C]=O(3394)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-765.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([355,410,600,1181,1341,1420,3056,375,486,586,623,710,938,1161,1294,1855,455,950,198.068,198.181,198.223],'cm^-1')),
        HinderedRotor(inertia=(0.265285,'amu*angstrom^2'), symmetry=1, barrier=(7.39092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265295,'amu*angstrom^2'), symmetry=1, barrier=(7.39114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973257,'amu*angstrom^2'), symmetry=1, barrier=(27.1365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.460061,0.0836366,-0.000133177,1.06851e-07,-3.36959e-11,-91997.3,27.8612], Tmin=(100,'K'), Tmax=(780.729,'K')), NASAPolynomial(coeffs=[12.9286,0.0197592,-1.0459e-05,2.06893e-09,-1.45551e-13,-93944.4,-29.207], Tmin=(780.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-765.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFOO) + group(Cds-OdCsH) + radical(O2sj(Cs-F1sF1sO2s)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]C(F)(F)C([O])(F)C=O(3269)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-703.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.926,'amu*angstrom^2'), symmetry=1, barrier=(21.2906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928319,'amu*angstrom^2'), symmetry=1, barrier=(21.3439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615782,0.0817456,-0.00012653,9.6194e-08,-2.71925e-11,-84439.1,27.3636], Tmin=(100,'K'), Tmax=(655.547,'K')), NASAPolynomial(coeffs=[12.0905,0.0216191,-1.15792e-05,2.30508e-09,-1.62854e-13,-86156,-24.7702], Tmin=(655.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-703.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = '[O]C(F)(F)OC=[C]F(3432)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7 C u0 p0 c0 {4,S} {8,D} {9,S}
8 C u1 p0 c0 {3,S} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-447.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05763,'amu*angstrom^2'), symmetry=1, barrier=(24.3169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06008,'amu*angstrom^2'), symmetry=1, barrier=(24.3732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612209,0.0686918,-8.45142e-05,4.79547e-08,-1.0322e-11,-53682.3,24.8238], Tmin=(100,'K'), Tmax=(1154.65,'K')), NASAPolynomial(coeffs=[18.9039,0.00532493,-2.19474e-06,4.25585e-10,-3.12435e-14,-57906.4,-66.0533], Tmin=(1154.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-447.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(O2sj(Cs-F1sF1sO2s)) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'FC1=COC(F)(F)OO1(3287)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {8,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-791.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13352,0.0528429,-3.23634e-05,-4.90756e-09,8.02324e-12,-95064.8,20.6454], Tmin=(100,'K'), Tmax=(963.692,'K')), NASAPolynomial(coeffs=[17.0725,0.0104459,-3.35625e-06,6.10893e-10,-4.56303e-14,-99240.2,-61.3859], Tmin=(963.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-791.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(124trioxene)"""),
)

species(
    label = '[O][C](F)C1OC(F)(F)O1(3433)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
9  C u1 p0 c0 {3,S} {6,S} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-741.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95337,0.0549218,-4.30306e-05,-3.85951e-08,6.59251e-11,-89059,24.3361], Tmin=(100,'K'), Tmax=(465.548,'K')), NASAPolynomial(coeffs=[7.07437,0.0285955,-1.51512e-05,3.02561e-09,-2.14659e-13,-89727.4,1.4885], Tmin=(465.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-741.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFOO) + ring(Cs-O2s-Cs(F)(F)-O2s) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = '[O]C(F)(F)O[C]=C=O(3400)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
7 C u1 p0 c0 {3,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
"""),
    E0 = (-328.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,1685,370,2120,512.5,787.5,184.774,184.882,185.369],'cm^-1')),
        HinderedRotor(inertia=(0.839833,'amu*angstrom^2'), symmetry=1, barrier=(20.0875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83147,'amu*angstrom^2'), symmetry=1, barrier=(20.0664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.027,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.846321,0.0639062,-8.31156e-05,4.92454e-08,-1.09823e-11,-39443.3,25.95], Tmin=(100,'K'), Tmax=(1119.79,'K')), NASAPolynomial(coeffs=[18.1051,0.00225608,-5.33019e-07,7.9913e-11,-5.80054e-15,-43308.6,-59.2666], Tmin=(1119.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsFFOO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(O2sj(Cs-F1sF1sO2s)) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)(F)O[C]=C(O)F(3434)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-668.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,486,586,623,710,938,1161,1294,293,496,537,1218,1685,370,206.219,531.04,531.131,531.135],'cm^-1')),
        HinderedRotor(inertia=(0.812336,'amu*angstrom^2'), symmetry=1, barrier=(18.6772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812334,'amu*angstrom^2'), symmetry=1, barrier=(18.6772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812278,'amu*angstrom^2'), symmetry=1, barrier=(18.6759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.244793,0.0796647,-0.000102294,6.10071e-08,-1.3892e-11,-80302.8,29.3243], Tmin=(100,'K'), Tmax=(1087.73,'K')), NASAPolynomial(coeffs=[19.6062,0.00846499,-4.1083e-06,8.28998e-10,-6.08642e-14,-84514.8,-65.7116], Tmin=(1087.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-668.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(O2sj(Cs-F1sF1sO2s)) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)=[C]OC(O)(F)F(3435)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-793.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611692,0.078787,-0.000108247,7.26419e-08,-1.91409e-11,-95362.3,28.5315], Tmin=(100,'K'), Tmax=(929.711,'K')), NASAPolynomial(coeffs=[14.8042,0.017726,-9.73248e-06,2.00109e-09,-1.45797e-13,-98001.3,-38.9047], Tmin=(929.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-793.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C(F)O[CH][C](F)OF(3436)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u1 p0 c0 {4,S} {8,S} {10,S}
8  C u1 p0 c0 {1,S} {5,S} {7,S}
9  C u0 p0 c0 {2,S} {4,S} {6,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-492.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,395,473,707,1436,482,664,788,1296,1923,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0261425,0.0965399,-0.000168447,1.44729e-07,-4.83298e-11,-59129,29.9851], Tmin=(100,'K'), Tmax=(802.049,'K')), NASAPolynomial(coeffs=[13.6307,0.0206081,-1.18102e-05,2.3726e-09,-1.6691e-13,-61068.1,-31.3203], Tmin=(802.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-492.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-CsOsHH) + group(CsCFHO) + group(COFOO) + radical(CCsJOC(O)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)(F)OC=[C]OF(3437)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {3,S} {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
8  C u0 p0 c0 {4,S} {9,D} {10,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-352.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,375,486,586,623,710,938,1161,1294,3010,987.5,1337.5,450,1655,1685,370,230.63,230.638,230.638,230.639],'cm^-1')),
        HinderedRotor(inertia=(0.586996,'amu*angstrom^2'), symmetry=1, barrier=(22.1596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58699,'amu*angstrom^2'), symmetry=1, barrier=(22.1597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587011,'amu*angstrom^2'), symmetry=1, barrier=(22.1596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0467288,0.0786164,-9.08786e-05,4.32028e-08,-6.07755e-12,-42258.5,30.4822], Tmin=(100,'K'), Tmax=(958.69,'K')), NASAPolynomial(coeffs=[22.7046,0.00236223,-1.74127e-07,1.95648e-11,-3.84793e-15,-47443,-82.2546], Tmin=(958.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-352.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(CsFFOO) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(O2sj(Cs-F1sF1sO2s)) + radical(C=CJO)"""),
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
    E0 = (-306.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (86.1621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (95.4461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-298.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-122.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-136.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-201.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-252.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-182.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-296.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (7.48444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-6.6218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (79.5018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-68.9612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (52.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-116.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-110.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (258.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-298.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-183.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (143.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-55.1439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-160.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (57.6648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (162.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['CF2O(48)', 'O=CC(=O)F(711)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C([O])(F)F(150)', '[CH]C(=O)F-2(3116)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', 'O=C(F)[CH]O[C](F)F(2827)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['O=C(F)C1OC(F)(F)O1(3354)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['[O]C(F)(F)OC1O[C]1F(3426)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['F[C]1[CH]OC(F)(F)OO1(3281)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.32511e+10,'s^-1'), n=0.274726, Ea=(169.057,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 167.9 to 169.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['[O]C1(F)[CH]OC(F)(F)O1(3308)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.16887e+10,'s^-1'), n=0.492453, Ea=(104.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_O] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C](F)F(237)', 'O=CC(=O)F(711)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(61.188,'m^3/(mol*s)'), n=1.29312, Ea=(22.4548,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.28976280499384915, var=2.1569028208455543, Tref=1000.0, N=37, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'O=C(F)[CH]OC(=O)F(3427)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(38.0644,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CF2O(48)', '[O][CH]C(=O)F(753)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.30416e-16,'m^3/(mol*s)'), n=6.08142, Ea=(72.862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.42494241855262277, var=2.081425725945606, Tref=1000.0, N=380, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]C(F)(F)OC=C=O(3428)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.9385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C](F)F(237)', '[O][CH]C(=O)F(753)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C(F)[CH]O[C](F)OF(3429)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(79.218,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C](F)OC(F)C(=O)F(3430)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(220.193,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C]=COC(F)(F)OF(3431)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(90.0713,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)(F)OC(F)[C]=O(3394)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(186.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3269)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(129.549,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(6)', '[O]C(F)(F)OC=[C]F(3432)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['FC1=COC(F)(F)OO1(3287)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['[O][C](F)C1OC(F)(F)O1(3433)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C(F)(F)O[C]=C=O(3400)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(290.474,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(F)(F)O[C]=C(O)F(3434)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    products = ['[O]C(F)=[C]OC(O)(F)F(3435)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C(F)O[CH][C](F)OF(3436)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(87.1002,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)(F)OC=[C]OF(3437)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3270)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(52.0458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1378',
    isomers = [
        '[O]C(F)(F)O[CH]C(=O)F(3270)',
    ],
    reactants = [
        ('CF2O(48)', 'O=CC(=O)F(711)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1378',
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

