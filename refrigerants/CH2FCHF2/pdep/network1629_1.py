species(
    label = '[O]OC(F)=C(F)OC(O)(F)F(5347)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1085.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,262,390,483,597,572,732,631,807,1275,1439,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.377738,'amu*angstrom^2'), symmetry=1, barrier=(8.68494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.455661,'amu*angstrom^2'), symmetry=1, barrier=(10.4765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378934,'amu*angstrom^2'), symmetry=1, barrier=(8.71243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.380318,'amu*angstrom^2'), symmetry=1, barrier=(8.74425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.763435,0.115187,-0.000202195,1.76455e-07,-5.95913e-11,-130413,34.4322], Tmin=(100,'K'), Tmax=(821.627,'K')), NASAPolynomial(coeffs=[14.2536,0.027419,-1.52004e-05,3.01343e-09,-2.10183e-13,-132386,-32.055], Tmin=(821.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1085.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = 'CF2O(49)',
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
    label = '[O]OC(F)C(=O)F(4704)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-582.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.519617,'amu*angstrom^2'), symmetry=1, barrier=(11.947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54035,'amu*angstrom^2'), symmetry=1, barrier=(35.4156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3788.12,'J/mol'), sigma=(5.7366,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=591.70 K, Pc=45.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66217,0.0558669,-8.74757e-05,7.36244e-08,-2.48646e-11,-69953.7,20.9883], Tmin=(100,'K'), Tmax=(757.994,'K')), NASAPolynomial(coeffs=[8.40288,0.0184184,-9.65371e-06,1.91163e-09,-1.34843e-13,-70921.6,-9.30828], Tmin=(757.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(C(=O)F)C(O)(F)F(8108)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {11,D}
8  O u1 p2 c0 {5,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
11 C u0 p0 c0 {4,S} {7,D} {9,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1237.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.824943,0.115019,-0.000195272,1.65173e-07,-5.45178e-11,-148636,33.5876], Tmin=(100,'K'), Tmax=(799.144,'K')), NASAPolynomial(coeffs=[15.571,0.0257457,-1.41792e-05,2.81715e-09,-1.9738e-13,-151027,-40.3976], Tmin=(799.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1237.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(COCsFO) + radical(ROOJ)"""),
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
    label = 'O=C(F)[C](F)OC(O)(F)F(8109)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
9  C u1 p0 c0 {3,S} {5,S} {10,S}
10 C u0 p0 c0 {4,S} {7,D} {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1305.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,280,501,1494,1531,611,648,830,1210,1753,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.473024,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223005,'amu*angstrom^2'), symmetry=1, barrier=(5.12733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282776,0.0893388,-0.000144291,1.2417e-07,-4.27218e-11,-156918,29.0586], Tmin=(100,'K'), Tmax=(762.103,'K')), NASAPolynomial(coeffs=[10.8483,0.0289362,-1.56648e-05,3.13229e-09,-2.21819e-13,-158385,-18.1007], Tmin=(762.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1305.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFOO) + group(COCsFO) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'OC(F)(F)OC1(F)OO[C]1F(8110)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {11,S}
8  O u0 p2 c0 {10,S} {12,S}
9  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1092.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.214575,0.101965,-0.000159214,1.24587e-07,-3.72714e-11,-131222,32.6883], Tmin=(100,'K'), Tmax=(658.581,'K')), NASAPolynomial(coeffs=[13.5601,0.0286362,-1.57342e-05,3.17136e-09,-2.25894e-13,-133261,-29.7148], Tmin=(658.581,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1092.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCFOO) + group(CsCFHO) + group(CsFFOO) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'OC(F)(F)O[C](F)C1(F)OO1(8111)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {10,S} {11,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {9,S}
8  O u0 p2 c0 {10,S} {12,S}
9  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
11 C u1 p0 c0 {4,S} {5,S} {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1102.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.378747,0.104477,-0.000166896,1.369e-07,-4.45529e-11,-132503,32.5856], Tmin=(100,'K'), Tmax=(754.568,'K')), NASAPolynomial(coeffs=[13.986,0.0283195,-1.54842e-05,3.11048e-09,-2.2102e-13,-134670,-32.6692], Tmin=(754.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1102.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCFOO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + ring(O2s-O2s-Cs(C-F)) + radical(CsCsF1sO2s)"""),
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
    label = '[O]OC(F)=C(F)O[C](O)F(8112)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,S} {11,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u0 p0 c0 {2,S} {5,S} {8,D}
10 C u1 p0 c0 {3,S} {4,S} {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-647.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,262,390,483,597,572,732,631,807,1275,1439,482,586,761,1411,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.204218,'amu*angstrom^2'), symmetry=1, barrier=(4.69537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204183,'amu*angstrom^2'), symmetry=1, barrier=(4.69456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204358,'amu*angstrom^2'), symmetry=1, barrier=(4.69859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704401,'amu*angstrom^2'), symmetry=1, barrier=(16.1956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.161891,0.101255,-0.000178707,1.57468e-07,-5.35275e-11,-77717.6,33.7979], Tmin=(100,'K'), Tmax=(830.647,'K')), NASAPolynomial(coeffs=[12.3354,0.0255505,-1.39662e-05,2.75053e-09,-1.90976e-13,-79258.2,-20.952], Tmin=(830.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-647.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFHOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]O[C]=C(F)OC(O)(F)F(8113)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
9  C u0 p0 c0 {3,S} {4,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-658.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,293,496,537,1218,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.174138,'amu*angstrom^2'), symmetry=1, barrier=(4.00377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174152,'amu*angstrom^2'), symmetry=1, barrier=(4.00409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174101,'amu*angstrom^2'), symmetry=1, barrier=(4.00292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60836,'amu*angstrom^2'), symmetry=1, barrier=(13.9874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.15858,0.10167,-0.000181555,1.61308e-07,-5.51227e-11,-79099.6,34.6705], Tmin=(100,'K'), Tmax=(833.927,'K')), NASAPolynomial(coeffs=[12.0728,0.0258135,-1.41935e-05,2.7978e-09,-1.94169e-13,-80541.9,-18.534], Tmin=(833.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]OC(F)=[C]OC(O)(F)F(8114)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-658.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,293,496,537,1218,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.174136,'amu*angstrom^2'), symmetry=1, barrier=(4.00373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174151,'amu*angstrom^2'), symmetry=1, barrier=(4.00408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1741,'amu*angstrom^2'), symmetry=1, barrier=(4.0029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608366,'amu*angstrom^2'), symmetry=1, barrier=(13.9875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.158677,0.101672,-0.00018156,1.61315e-07,-5.51265e-11,-79099.6,34.6708], Tmin=(100,'K'), Tmax=(833.901,'K')), NASAPolynomial(coeffs=[12.073,0.0258132,-1.41933e-05,2.79775e-09,-1.94164e-13,-80542,-18.535], Tmin=(833.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'O[C](F)F(2562)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-470.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,493,600,700,1144,1293],'cm^-1')),
        HinderedRotor(inertia=(0.307072,'amu*angstrom^2'), symmetry=1, barrier=(7.06019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94331,0.003562,5.16413e-05,-1.18327e-07,7.89485e-11,-56559.5,10.5484], Tmin=(10,'K'), Tmax=(519.044,'K')), NASAPolynomial(coeffs=[4.2538,0.012858,-9.00307e-06,2.95251e-09,-3.63989e-13,-56749.2,7.73732], Tmin=(519.044,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-470.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]O[C](F)C(=O)F(4794)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u1 p0 c0 {1,S} {3,S} {7,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
"""),
    E0 = (-443.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,280,501,1494,1531,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(4.11501,'amu*angstrom^2'), symmetry=1, barrier=(94.6122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.22048,'amu*angstrom^2'), symmetry=1, barrier=(97.0371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.016,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93926,0.0519595,-9.19904e-05,8.61301e-08,-3.09897e-11,-53220.8,20.2084], Tmin=(100,'K'), Tmax=(837.273,'K')), NASAPolynomial(coeffs=[5.80387,0.0200068,-1.05786e-05,2.06414e-09,-1.42885e-13,-53395.1,5.07392], Tmin=(837.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(O)(F)F(2944)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-638.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,486,586,623,710,938,1161,1294],'cm^-1')),
        HinderedRotor(inertia=(0.596909,'amu*angstrom^2'), symmetry=1, barrier=(13.7241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90705,0.00545186,7.4862e-05,-1.61389e-07,9.87578e-11,-76836.8,10.4258], Tmin=(10,'K'), Tmax=(579.229,'K')), NASAPolynomial(coeffs=[5.15887,0.0182822,-1.39769e-05,4.86853e-09,-6.26572e-13,-77342.1,1.96031], Tmin=(579.229,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-638.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]C(O)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)=[C]F(4600)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u1 p0 c0 {2,S} {5,D}
"""),
    E0 = (0.671019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,293,496,537,1218,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.147304,'amu*angstrom^2'), symmetry=1, barrier=(3.38682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94084,0.0545835,-0.000112998,1.10274e-07,-3.96452e-11,145.962,18.1125], Tmin=(100,'K'), Tmax=(865.493,'K')), NASAPolynomial(coeffs=[5.76978,0.015686,-8.83972e-06,1.73986e-09,-1.19296e-13,277.257,4.78062], Tmin=(865.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.671019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
)

species(
    label = 'OH(7)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]OC(F)=C(F)O[C](F)F(8115)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {1,S} {6,S} {8,D}
10 C u1 p0 c0 {3,S} {4,S} {5,S}
"""),
    E0 = (-658.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,493,600,700,1144,1293,180,180,1752.52],'cm^-1')),
        HinderedRotor(inertia=(0.346149,'amu*angstrom^2'), symmetry=1, barrier=(7.95864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345809,'amu*angstrom^2'), symmetry=1, barrier=(7.95083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346117,'amu*angstrom^2'), symmetry=1, barrier=(7.95792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.253901,0.0920227,-0.000164534,1.48463e-07,-5.18485e-11,-79028.1,32.5253], Tmin=(100,'K'), Tmax=(816.236,'K')), NASAPolynomial(coeffs=[10.6317,0.0257956,-1.45824e-05,2.9192e-09,-2.04969e-13,-80210.2,-12.2981], Tmin=(816.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-Cd))"""),
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
    label = '[O]OC(F)=C(F)OC([O])(F)F(8116)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u1 p2 c0 {9,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {3,S} {6,S} {10,D}
"""),
    E0 = (-819.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,486,586,623,710,938,1161,1294,262,390,483,597,572,732,631,807,1275,1439,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.303397,'amu*angstrom^2'), symmetry=1, barrier=(6.97569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303704,'amu*angstrom^2'), symmetry=1, barrier=(6.98275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303886,'amu*angstrom^2'), symmetry=1, barrier=(6.98694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (176.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.142849,0.101358,-0.000160407,1.18664e-07,-3.03178e-11,-98380.9,33.1473], Tmin=(100,'K'), Tmax=(624.813,'K')), NASAPolynomial(coeffs=[14.3564,0.0243929,-1.37076e-05,2.75919e-09,-1.95402e-13,-100502,-32.4614], Tmin=(624.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-819.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-F1sF1sO2s)) + radical(ROOJ)"""),
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
    label = 'OC(F)(F)OC(F)=[C]F(4072)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-890.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,293,496,537,1218,167,640,1190,316.868,317.799,318.993,319.38],'cm^-1')),
        HinderedRotor(inertia=(0.0874035,'amu*angstrom^2'), symmetry=1, barrier=(6.26193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201309,'amu*angstrom^2'), symmetry=1, barrier=(14.5019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.497026,'amu*angstrom^2'), symmetry=1, barrier=(35.3422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446696,0.0830247,-0.000120076,8.48023e-08,-2.34473e-11,-106967,25.9266], Tmin=(100,'K'), Tmax=(888.49,'K')), NASAPolynomial(coeffs=[14.9691,0.0176444,-9.69769e-06,1.9812e-09,-1.43375e-13,-109548,-42.4188], Tmin=(888.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(Cdj(Cd-F1sO2s)(F1s))"""),
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
    label = '[O]OC(F)=C(F)OC(=O)F(8117)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u0 p0 c0 {1,S} {5,S} {8,D}
10 C u0 p0 c0 {3,S} {4,S} {6,D}
"""),
    E0 = (-794.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,262,390,483,597,572,732,631,807,1275,1439,482,664,788,1296,1923,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.23214,'amu*angstrom^2'), symmetry=1, barrier=(5.33736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23211,'amu*angstrom^2'), symmetry=1, barrier=(5.33667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231723,'amu*angstrom^2'), symmetry=1, barrier=(5.32777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15231,0.0990773,-0.000193657,1.85803e-07,-6.71465e-11,-95442,31.0166], Tmin=(100,'K'), Tmax=(841.783,'K')), NASAPolynomial(coeffs=[8.07967,0.0309218,-1.7885e-05,3.5734e-09,-2.49032e-13,-95696.5,0.55269], Tmin=(841.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-794.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(COFOO) + radical(ROOJ)"""),
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
    label = '[O]OC#COC(O)(F)F(8118)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
8  C u0 p0 c0 {3,S} {9,T}
9  C u0 p0 c0 {5,S} {8,T}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-471.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,2100,2250,500,550,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11179,'amu*angstrom^2'), symmetry=1, barrier=(25.5623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11181,'amu*angstrom^2'), symmetry=1, barrier=(25.5626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11204,'amu*angstrom^2'), symmetry=1, barrier=(25.568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11077,'amu*angstrom^2'), symmetry=1, barrier=(25.5387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0567692,0.0856668,-0.00011757,7.47482e-08,-1.81208e-11,-56556.8,27.8219], Tmin=(100,'K'), Tmax=(1023.31,'K')), NASAPolynomial(coeffs=[19.8604,0.00825555,-4.09585e-06,8.21108e-10,-5.97159e-14,-60609.8,-68.1751], Tmin=(1023.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-471.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsCt) + group(O2s-OsH) + group(CsFFOO) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=C=O(2248)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u1 p2 c0 {2,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (-105.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,197,221,431,657,2120,512.5,787.5,2787.33],'cm^-1')),
        HinderedRotor(inertia=(0.0346032,'amu*angstrom^2'), symmetry=1, barrier=(3.05444,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37136,0.043506,-8.64301e-05,8.52867e-08,-3.11222e-11,-12657.6,18.7731], Tmin=(100,'K'), Tmax=(866.808,'K')), NASAPolynomial(coeffs=[4.28985,0.0166729,-8.88159e-06,1.71369e-09,-1.1673e-13,-12314.7,13.688], Tmin=(866.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = '[O]C(F)(F)OC(F)=C(F)OO(8119)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {7,S} {11,S}
7  O u0 p2 c0 {6,S} {12,S}
8  O u1 p2 c0 {9,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {3,S} {6,S} {10,D}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-971.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,375,486,586,623,710,938,1161,1294,262,390,483,597,572,732,631,807,1275,1439,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.214081,'amu*angstrom^2'), symmetry=1, barrier=(4.92215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94392,'amu*angstrom^2'), symmetry=1, barrier=(21.7026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214153,'amu*angstrom^2'), symmetry=1, barrier=(4.9238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96841,'amu*angstrom^2'), symmetry=1, barrier=(45.2576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42711,0.105631,-0.000161054,1.23069e-07,-3.7202e-11,-116650,33.3861], Tmin=(100,'K'), Tmax=(811.139,'K')), NASAPolynomial(coeffs=[15.5805,0.0266916,-1.50743e-05,3.08993e-09,-2.23074e-13,-119247,-40.4909], Tmin=(811.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-971.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(O2sj(Cs-F1sF1sO2s))"""),
)

species(
    label = 'O[C](F)OC(F)=C(F)OOF(8120)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,S} {12,S}
8  O u0 p2 c0 {4,S} {6,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {1,S} {6,S} {9,D}
11 C u1 p0 c0 {3,S} {5,S} {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-704.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,262,390,483,597,572,732,631,807,1275,1439,482,586,761,1411],'cm^-1')),
        HinderedRotor(inertia=(0.70588,'amu*angstrom^2'), symmetry=1, barrier=(16.2296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708568,'amu*angstrom^2'), symmetry=1, barrier=(16.2914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86246,'amu*angstrom^2'), symmetry=1, barrier=(42.8216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85946,'amu*angstrom^2'), symmetry=1, barrier=(42.7527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710118,'amu*angstrom^2'), symmetry=1, barrier=(16.327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744803,0.11403,-0.000195126,1.68573e-07,-5.71683e-11,-84516.4,36.5134], Tmin=(100,'K'), Tmax=(781.904,'K')), NASAPolynomial(coeffs=[14.4282,0.0285968,-1.6244e-05,3.27658e-09,-2.31953e-13,-86650.4,-31.4278], Tmin=(781.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2sFO) + group(CsFHOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsF1sO2sO2s)"""),
)

species(
    label = 'OC(F)(F)O[C]=C(F)OOF(8121)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {4,S} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u1 p0 c0 {5,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-715.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3615,1310,387.5,850,1000,277,555,632,435,565,619,662,854,1178,1396,293,496,537,1218,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.74075,'amu*angstrom^2'), symmetry=1, barrier=(40.0233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667954,'amu*angstrom^2'), symmetry=1, barrier=(15.3576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66779,'amu*angstrom^2'), symmetry=1, barrier=(15.3538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74083,'amu*angstrom^2'), symmetry=1, barrier=(40.0252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668165,'amu*angstrom^2'), symmetry=1, barrier=(15.3624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.748921,0.114542,-0.000198348,1.72951e-07,-5.90148e-11,-85898.1,37.412], Tmin=(100,'K'), Tmax=(790.445,'K')), NASAPolynomial(coeffs=[14.1889,0.0288178,-1.6446e-05,3.31767e-09,-2.34621e-13,-87943.1,-29.1392], Tmin=(790.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2sFO) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=CJO)"""),
)

species(
    label = 'OC(F)(F)OC(F)=[C]OOF(8122)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u0 p2 c0 {4,S} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
10 C u0 p0 c0 {3,S} {5,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-715.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3615,1310,387.5,850,1000,277,555,632,435,565,619,662,854,1178,1396,293,496,537,1218,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.74166,'amu*angstrom^2'), symmetry=1, barrier=(40.0442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66774,'amu*angstrom^2'), symmetry=1, barrier=(15.3527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667324,'amu*angstrom^2'), symmetry=1, barrier=(15.3431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73944,'amu*angstrom^2'), symmetry=1, barrier=(39.9932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668882,'amu*angstrom^2'), symmetry=1, barrier=(15.3789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.749454,0.114549,-0.000198378,1.72998e-07,-5.90391e-11,-85898.1,37.4139], Tmin=(100,'K'), Tmax=(790.224,'K')), NASAPolynomial(coeffs=[14.1897,0.0288163,-1.6445e-05,3.31744e-09,-2.34601e-13,-87943.4,-29.1439], Tmin=(790.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2sFO) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
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
    E0 = (-460.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-467.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-540.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-437.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-517.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-52.1765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-63.6623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-63.6623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-391.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-115.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-107.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-85.0599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-376.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-389.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (42.0466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-426.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-376.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-122.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-130.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-108.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    products = ['CF2O(49)', '[O]OC(F)C(=O)F(4704)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(103.065,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    products = ['[O]OC(F)(C(=O)F)C(O)(F)F(8108)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(95.5161,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', 'O=C(F)[C](F)OC(O)(F)F(8109)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    products = ['OC(F)(F)OC1(F)OO[C]1F(8110)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    products = ['OC(F)(F)O[C](F)C1(F)OO1(8111)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]OC(F)=C(F)O[C](O)F(8112)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', '[O]O[C]=C(F)OC(O)(F)F(8113)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', '[O]OC(F)=[C]OC(O)(F)F(8114)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O[C](F)F(2562)', '[O]O[C](F)C(=O)F(4794)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(O)(F)F(2944)', '[O]OC(F)=[C]F(4600)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OH(7)', '[O]OC(F)=C(F)O[C](F)F(8115)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[O]OC(F)=C(F)OC([O])(F)F(8116)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O2(2)', 'OC(F)(F)OC(F)=[C]F(4072)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', '[O]OC(F)=C(F)OC(=O)F(8117)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(163.71,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F2(78)', '[O]OC#COC(O)(F)F(8118)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    products = ['HF(38)', 'CF2O(49)', '[O]OC(F)=C=O(2248)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.14547e+10,'s^-1'), n=0.631481, Ea=(137,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)(F)OC(F)=C(F)OO(8119)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(627096,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;O_rad_out;XH_out] for rate rule [R7H;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O[C](F)OC(F)=C(F)OOF(8120)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(59.5649,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OC(F)(F)O[C]=C(F)OOF(8121)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(62.3079,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OC(F)(F)OC(F)=[C]OOF(8122)'],
    products = ['[O]OC(F)=C(F)OC(O)(F)F(5347)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(84.7403,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1629',
    isomers = [
        '[O]OC(F)=C(F)OC(O)(F)F(5347)',
    ],
    reactants = [
        ('CF2O(49)', '[O]OC(F)C(=O)F(4704)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1629',
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

