species(
    label = 'OC(F)(F)C=C=CF(9336)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {3,S} {8,D} {10,S}
8  C u0 p0 c0 {6,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-647.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,113,247,382,1207,3490,540,610,2055,223.608,1057.29],'cm^-1')),
        HinderedRotor(inertia=(0.128562,'amu*angstrom^2'), symmetry=1, barrier=(4.54652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130093,'amu*angstrom^2'), symmetry=1, barrier=(4.53955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02751,0.0706022,-9.70746e-05,7.188e-08,-2.15372e-11,-77749.6,23.2974], Tmin=(100,'K'), Tmax=(812.25,'K')), NASAPolynomial(coeffs=[10.3252,0.0248146,-1.25176e-05,2.47842e-09,-1.7628e-13,-79260.1,-19.6252], Tmin=(812.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-647.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds)"""),
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
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'OC=C=CF(11416)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u0 p0 c0 {2,S} {5,D} {6,S}
4 C u0 p0 c0 {1,S} {5,D} {7,S}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-160.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,113,247,382,1207,3490,540,610,2055,581.08],'cm^-1')),
        HinderedRotor(inertia=(0.136464,'amu*angstrom^2'), symmetry=1, barrier=(32.6977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90003,0.00602671,9.90309e-05,-2.10076e-07,1.30868e-10,-19266.1,9.79642], Tmin=(10,'K'), Tmax=(549.647,'K')), NASAPolynomial(coeffs=[4.01944,0.0272351,-1.90972e-05,6.27839e-09,-7.75988e-13,-19612.7,6.25808], Tmin=(549.647,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-160.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""OCDCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OF(169)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]C=C=CF-2(11687)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {5,D} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {5,D} {8,S}
5 C u0 p0 c0 {3,D} {4,D}
6 C u0 p1 c0 {2,S} {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (135.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,113,247,382,1207,3490,540,610,2055,315,622,1128,547.186],'cm^-1')),
        HinderedRotor(inertia=(0.000563248,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44619,0.035976,-3.18673e-05,1.43396e-08,-2.67643e-12,16395.3,16.7517], Tmin=(100,'K'), Tmax=(1239.94,'K')), NASAPolynomial(coeffs=[8.49186,0.016473,-8.27398e-06,1.65443e-09,-1.18832e-13,14896,-13.7155], Tmin=(1239.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + group(CJ2_singlet-FC)"""),
)

species(
    label = 'OC(F)=CC(F)=CF(11688)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {11,S}
5  C u0 p0 c0 {6,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {2,S} {4,S} {5,D}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-662.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.692856,0.0710931,-8.24183e-05,4.66973e-08,-1.03222e-11,-79504.3,22.4786], Tmin=(100,'K'), Tmax=(1108.1,'K')), NASAPolynomial(coeffs=[16.1056,0.0154557,-7.10252e-06,1.38433e-09,-9.8963e-14,-82920,-53.4606], Tmin=(1108.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-662.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'OC(C=C(F)F)=CF(11689)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {3,S} {5,D} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-677.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674331,0.0714025,-8.51638e-05,4.99845e-08,-1.14116e-11,-81398.1,22.5002], Tmin=(100,'K'), Tmax=(1076.64,'K')), NASAPolynomial(coeffs=[15.8935,0.014859,-6.38571e-06,1.20394e-09,-8.45344e-14,-84675.2,-52.0471], Tmin=(1076.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-677.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFF)"""),
)

species(
    label = 'OC(F)(F)[C]=C[CH]F(11690)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {7,S} {8,D} {9,S}
7  C u1 p0 c0 {3,S} {6,S} {10,S}
8  C u1 p0 c0 {5,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-456.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,140,424,499,621,667,843,876,1082,3010,987.5,1337.5,450,1655,234,589,736,816,1240,3237,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.364375,'amu*angstrom^2'), symmetry=1, barrier=(8.3777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36533,'amu*angstrom^2'), symmetry=1, barrier=(8.39967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992331,'amu*angstrom^2'), symmetry=1, barrier=(22.8157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71696,0.0786261,-0.000123239,1.04526e-07,-3.56262e-11,-54797.5,25.4295], Tmin=(100,'K'), Tmax=(757.122,'K')), NASAPolynomial(coeffs=[9.95078,0.0266614,-1.39853e-05,2.77585e-09,-1.9614e-13,-56104.6,-15.9471], Tmin=(757.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'OC(F)(F)[CH]C=[C]F(11691)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u1 p0 c0 {5,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-461.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05777,'amu*angstrom^2'), symmetry=1, barrier=(47.3121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934425,'amu*angstrom^2'), symmetry=1, barrier=(21.4843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.933344,'amu*angstrom^2'), symmetry=1, barrier=(21.4594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707681,0.0755386,-9.56833e-05,6.1067e-08,-1.54463e-11,-55392.2,23.5325], Tmin=(100,'K'), Tmax=(964.329,'K')), NASAPolynomial(coeffs=[13.9383,0.020659,-1.03196e-05,2.05341e-09,-1.47352e-13,-57944,-39.8172], Tmin=(964.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-461.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + radical(C=CCJCO) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C(F)(F)C[C]=CF(8611)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-321.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,351,323,533,609,664,892,1120,1201,615,860,1140,1343,3152,1685,370,180,180,3413.54],'cm^-1')),
        HinderedRotor(inertia=(0.299675,'amu*angstrom^2'), symmetry=1, barrier=(6.89011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960365,'amu*angstrom^2'), symmetry=1, barrier=(22.0807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.56,'J/mol'), sigma=(6.00584,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.90 K, Pc=37.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32258,0.0628215,-7.30405e-05,4.52866e-08,-1.14764e-11,-38617.1,26.1899], Tmin=(100,'K'), Tmax=(948.796,'K')), NASAPolynomial(coeffs=[10.4845,0.0241958,-1.19751e-05,2.3792e-09,-1.70605e-13,-40355.6,-17.5295], Tmin=(948.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + radical(O2sj(Cs-CsF1sF1s)) + radical(Cdj(Cs-CsHH)(Cd-F1sH))"""),
)

species(
    label = '[O]C(F)(F)[CH]C=CF(9340)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u1 p0 c0 {5,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u0 p0 c0 {3,S} {7,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-451.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,280.654,281.363,281.56],'cm^-1')),
        HinderedRotor(inertia=(0.456867,'amu*angstrom^2'), symmetry=1, barrier=(25.8016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739558,'amu*angstrom^2'), symmetry=1, barrier=(41.6658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.95718,0.0685522,-7.71485e-05,4.39193e-08,-9.99143e-12,-54215.2,22.7311], Tmin=(100,'K'), Tmax=(1062.98,'K')), NASAPolynomial(coeffs=[13.5235,0.0212646,-1.04195e-05,2.06881e-09,-1.48643e-13,-56886.8,-38.6618], Tmin=(1062.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-451.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + radical(O2sj(Cs-CsF1sF1s)) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C(F)(F)C=[C]CF(11124)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u1 p0 c0 {6,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-339.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,122,584,932,1043,1247,1392,1461,2986,3039,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.788029,'amu*angstrom^2'), symmetry=1, barrier=(18.1183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78801,'amu*angstrom^2'), symmetry=1, barrier=(18.1179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.56,'J/mol'), sigma=(6.00584,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.90 K, Pc=37.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49952,0.0562652,-5.25294e-05,2.48156e-08,-4.78682e-12,-40742.3,26.1114], Tmin=(100,'K'), Tmax=(1220.73,'K')), NASAPolynomial(coeffs=[11.7753,0.022594,-1.11549e-05,2.21992e-09,-1.59286e-13,-43251,-25.5127], Tmin=(1220.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(O2sj(Cs-F1sF1sCd)) + radical(Cdj(Cs-F1sHH)(Cd-CsH))"""),
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
    label = 'OC(F)=C[C]=CF(11692)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {10,S}
4  C u0 p0 c0 {5,D} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {4,D}
6  C u0 p0 c0 {2,S} {7,D} {9,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-268.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,615,860,1140,1343,3152,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.944083,'amu*angstrom^2'), symmetry=1, barrier=(21.7063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939778,'amu*angstrom^2'), symmetry=1, barrier=(21.6073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938084,0.0637582,-7.74294e-05,4.61208e-08,-1.05662e-11,-32135.7,22.3335], Tmin=(100,'K'), Tmax=(1081.07,'K')), NASAPolynomial(coeffs=[15.3174,0.0105546,-3.60928e-06,5.98355e-10,-3.91066e-14,-35244.7,-48.1594], Tmin=(1081.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFO) + group(CdCFH) + radical(Cdj(Cd-CdH)(Cd-F1sH))"""),
)

species(
    label = '[CH]=C=CC(O)(F)F(10716)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {8,S}
6  C u0 p0 c0 {5,D} {7,D}
7  C u1 p0 c0 {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-319.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.804984,'amu*angstrom^2'), symmetry=1, barrier=(18.5082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.804421,'amu*angstrom^2'), symmetry=1, barrier=(18.4952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18118,0.0645855,-8.90997e-05,6.40601e-08,-1.81821e-11,-38275.4,21.8556], Tmin=(100,'K'), Tmax=(865.567,'K')), NASAPolynomial(coeffs=[11.3359,0.0176629,-7.79291e-06,1.44355e-09,-9.86152e-14,-40033.5,-25.6702], Tmin=(865.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC=[C]C=C(F)F(11645)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 C u0 p0 c0 {3,S} {7,D} {9,S}
7 C u1 p0 c0 {4,S} {6,D}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-264.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,615,860,1140,1343,3152,1685,370,180,418.578],'cm^-1')),
        HinderedRotor(inertia=(1.32222,'amu*angstrom^2'), symmetry=1, barrier=(30.4004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61467,0.0369259,-2.62925e-05,6.40149e-09,3.99651e-13,-31763.4,12.4832], Tmin=(10,'K'), Tmax=(1021.64,'K')), NASAPolynomial(coeffs=[9.31624,0.0209883,-1.22678e-05,3.3676e-09,-3.55e-13,-33261.6,-16.7766], Tmin=(1021.64,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-264.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FCD[C]CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C=C=CF(9337)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {3,S} {8,D} {10,S}
8  C u0 p0 c0 {6,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-399.285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,3010,987.5,1337.5,450,1655,113,247,382,1207,3490,540,610,2055,183.86,184.849,1214.98],'cm^-1')),
        HinderedRotor(inertia=(0.982072,'amu*angstrom^2'), symmetry=1, barrier=(23.632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40988,0.0573077,-6.02593e-05,3.18937e-08,-6.75481e-12,-47929.8,23.4255], Tmin=(100,'K'), Tmax=(1136.58,'K')), NASAPolynomial(coeffs=[12.4715,0.018378,-8.88157e-06,1.75771e-09,-1.26122e-13,-50444.3,-31.3566], Tmin=(1136.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(O2sj(Cs-F1sF1sCd))"""),
)

species(
    label = 'O[C](F)F(4999)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94331,0.003562,5.16413e-05,-1.18327e-07,7.89485e-11,-56559.5,10.5484], Tmin=(10,'K'), Tmax=(519.044,'K')), NASAPolynomial(coeffs=[4.2538,0.012858,-9.00307e-06,2.95251e-09,-3.63989e-13,-56749.2,7.73732], Tmin=(519.044,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-470.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C#C[CH]F(4374)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {5,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {6,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (160.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([300,424,1135,1289,3214,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.63779,'amu*angstrom^2'), symmetry=1, barrier=(37.6559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0463,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90976,0.00602051,7.44069e-05,-1.92463e-07,1.43212e-10,19366.8,8.53447], Tmin=(10,'K'), Tmax=(475.298,'K')), NASAPolynomial(coeffs=[4.8953,0.0152178,-9.82007e-06,3.10302e-09,-3.78569e-13,19075.6,2.43445], Tmin=(475.298,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(160.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C#C[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)C#C[CH]F(11693)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u1 p0 c0 {3,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,T}
8  C u0 p0 c0 {6,S} {7,T}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-514.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,157,399,436,611,601,744,1127,300,424,1135,1289,3214,2100,2250,500,550,180,1784.61],'cm^-1')),
        HinderedRotor(inertia=(2.90925,'amu*angstrom^2'), symmetry=1, barrier=(66.8895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540486,'amu*angstrom^2'), symmetry=1, barrier=(12.4268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.90353,'amu*angstrom^2'), symmetry=1, barrier=(66.7578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.67018,0.0780362,-0.000124981,1.00804e-07,-3.12154e-11,-61789.1,24.8532], Tmin=(100,'K'), Tmax=(897.232,'K')), NASAPolynomial(coeffs=[11.9945,0.0176873,-7.60029e-06,1.33446e-09,-8.60944e-14,-63424.2,-26.3399], Tmin=(897.232,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-514.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CsCFHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Csj(F1s)(Ct-Ct)(H))"""),
)

species(
    label = 'OC(F)(F)C=C=[C]F(10118)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {9,S}
7  C u0 p0 c0 {6,D} {8,D}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-439.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,540,610,2055,137,207,812,180],'cm^-1')),
        HinderedRotor(inertia=(0.161232,'amu*angstrom^2'), symmetry=1, barrier=(3.70703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.389182,'amu*angstrom^2'), symmetry=1, barrier=(8.94805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12146,0.0704122,-0.000100982,6.89914e-08,-1.49515e-11,-52768.5,23.6361], Tmin=(100,'K'), Tmax=(615.452,'K')), NASAPolynomial(coeffs=[10.3182,0.0221339,-1.13302e-05,2.22346e-09,-1.5609e-13,-54118.2,-18.0372], Tmin=(615.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-439.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(Cdj(Cdd-Cd)(F1s))"""),
)

species(
    label = 'OC(F)(F)[C]C=CF(11694)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {7,D} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {6,D} {10,S}
8  C u0 p1 c0 {5,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-401.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.73529,'amu*angstrom^2'), symmetry=1, barrier=(39.8977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73755,'amu*angstrom^2'), symmetry=1, barrier=(39.9496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73229,'amu*angstrom^2'), symmetry=1, barrier=(39.8287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.535961,0.0839076,-0.000136876,1.19694e-07,-4.11796e-11,-48115.1,35.0759], Tmin=(100,'K'), Tmax=(818.575,'K')), NASAPolynomial(coeffs=[9.42161,0.0288641,-1.4712e-05,2.85401e-09,-1.97905e-13,-49180.4,-3.63472], Tmin=(818.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(F)=CC(O)(F)F-2(11695)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {9,S}
7  C u0 p0 c0 {3,S} {6,D} {8,S}
8  C u0 p1 c0 {7,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-401.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,405,756,886,1212,180,180,180,671.599,2646.68],'cm^-1')),
        HinderedRotor(inertia=(0.21609,'amu*angstrom^2'), symmetry=1, barrier=(4.96834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210116,'amu*angstrom^2'), symmetry=1, barrier=(4.83098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218658,'amu*angstrom^2'), symmetry=1, barrier=(5.02738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598065,0.0821121,-0.000133984,1.16969e-07,-4.07164e-11,-48194.4,24.3624], Tmin=(100,'K'), Tmax=(771.219,'K')), NASAPolynomial(coeffs=[9.90989,0.0273086,-1.47372e-05,2.94811e-09,-2.08767e-13,-49437.2,-16.8882], Tmin=(771.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCCF) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'OC(F)(F)C=C[C]F(11696)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {11,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {9,S}
7  C u0 p0 c0 {6,D} {8,S} {10,S}
8  C u0 p1 c0 {3,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-487.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,2995,3025,975,1000,1300,1375,400,500,1630,1680,315,622,1128,180],'cm^-1')),
        HinderedRotor(inertia=(0.202591,'amu*angstrom^2'), symmetry=1, barrier=(4.65797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.448603,'amu*angstrom^2'), symmetry=1, barrier=(10.3143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885568,'amu*angstrom^2'), symmetry=1, barrier=(20.3609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19765,0.0669595,-8.25359e-05,5.37137e-08,-1.42755e-11,-58503.1,23.769], Tmin=(100,'K'), Tmax=(906.092,'K')), NASAPolynomial(coeffs=[10.6911,0.0250522,-1.31633e-05,2.67444e-09,-1.93875e-13,-60223.5,-21.0954], Tmin=(906.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-487.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CJ2_singlet-FC)"""),
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
    label = 'O=C(F)C=C=CF(11697)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {5,S} {7,D} {8,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 C u0 p0 c0 {2,S} {7,D} {9,S}
7 C u0 p0 c0 {4,D} {6,D}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-353.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,255,533,799,832,1228,113,247,382,1207,3490,540,610,2055,1334.02,1337.72],'cm^-1')),
        HinderedRotor(inertia=(0.311791,'amu*angstrom^2'), symmetry=1, barrier=(7.16869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99533,0.0480688,-5.82148e-05,4.00132e-08,-1.15101e-11,-42415.4,19.3901], Tmin=(100,'K'), Tmax=(833.126,'K')), NASAPolynomial(coeffs=[7.4674,0.0217983,-1.09198e-05,2.17066e-09,-1.55314e-13,-43327.2,-6.01079], Tmin=(833.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(COCFO) + group(CdCddFH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'OC(F)=C=C=CF(11698)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {9,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u0 p0 c0 {2,S} {7,D} {8,S}
6 C u0 p0 c0 {4,D} {7,D}
7 C u0 p0 c0 {5,D} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-225.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,197,221,431,657,113,247,382,1207,3490,540,563.333,586.667,610,1970,2140,418.62,418.867],'cm^-1')),
        HinderedRotor(inertia=(0.000959343,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31374,0.0594792,-7.53264e-05,4.69489e-08,-1.14061e-11,-26986.2,19.7987], Tmin=(100,'K'), Tmax=(1011.56,'K')), NASAPolynomial(coeffs=[13.173,0.0125838,-5.78658e-06,1.11842e-09,-7.92921e-14,-29385.5,-37.5516], Tmin=(1011.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CdCddFO) + group(CdCddFH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    E0 = (-228.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (102.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (309.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-88.1892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-174.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-176.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-181.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (14.0185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-169.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-56.8896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (64.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (11.4538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (21.9199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (70.1465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-50.6456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-40.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (29.8832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-103.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-27.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-143.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-214.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-49.3319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)(F)C=C=CF(9336)'],
    products = ['CF2O(49)', 'C#CCF(4379)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(161.034,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'OC=C=CF(11416)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.50469e+53,'m^3/(mol*s)'), n=-13.541, Ea=(209.104,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction3',
    reactants = ['OF(169)', 'F[C]C=C=CF-2(11687)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(11.5949,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OC(F)(F)C=C=CF(9336)'],
    products = ['OC(F)=CC(F)=CF(11688)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.45932e+11,'s^-1'), n=0.63878, Ea=(301.483,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='F_Ext-2R!H-R',), comment="""Estimated from node F_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OC(F)(F)C=C=CF(9336)'],
    products = ['OC(C=C(F)F)=CF(11689)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(215.426,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OC(F)(F)[C]=C[CH]F(11690)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC(F)(F)[CH]C=[C]F(11691)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)(F)C[C]=CF(8611)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)(F)[CH]C=CF(9340)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)(F)C=[C]CF(11124)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'OC(F)=C[C]=CF(11692)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(2.18769,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[CH]=C=CC(O)(F)F(10716)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.11549e+09,'m^3/(mol*s)'), n=-0.68237, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.29431919638206655, var=1.0853977775937997, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_Sp-3BrBrCCOO=1C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_Sp-3BrBrCCOO=1C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OH(7)', 'FC=[C]C=C(F)F(11645)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', '[O]C(F)(F)C=C=CF(9337)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O[C](F)F(4999)', 'C#C[CH]F(4374)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(1.03291,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'OC(F)(F)C#C[CH]F(11693)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(4.56093,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'OC(F)(F)C=C=[C]F(10118)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -1.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['OC(F)(F)[C]C=CF(11694)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.33333e+12,'s^-1'), n=8.2394e-08, Ea=(40.013,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C(F)=CC(O)(F)F-2(11695)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(116.056,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OC(F)(F)C=C[C]F(11696)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.33334e+12,'s^-1'), n=-1.40567e-07, Ea=(86.5393,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_4CClFINOPSSi->F',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_4CClFINOPSSi->F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'O=C(F)C=C=CF(11697)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(162.546,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'OC(F)=C=C=CF(11698)'],
    products = ['OC(F)(F)C=C=CF(9336)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(199.334,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #4055',
    isomers = [
        'OC(F)(F)C=C=CF(9336)',
    ],
    reactants = [
        ('CF2O(49)', 'C#CCF(4379)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4055',
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

