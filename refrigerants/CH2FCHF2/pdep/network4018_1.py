species(
    label = 'F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {13,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
14 C u1 p0 c0 {7,S} {13,S} {16,S}
15 H u0 p0 c0 {13,S}
16 H u0 p0 c0 {14,S}
"""),
    E0 = (-1463.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (231.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84757,0.136448,-0.000184281,1.2245e-07,-3.20959e-11,-175757,42.7433], Tmin=(100,'K'), Tmax=(932.486,'K')), NASAPolynomial(coeffs=[22.2126,0.033239,-1.82592e-05,3.75525e-09,-2.73706e-13,-180244,-71.6512], Tmin=(932.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1463.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[O]OC(F)(F)C(F)(F)F(13321)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {9,S}
4 F u0 p3 c0 {9,S}
5 F u0 p3 c0 {9,S}
6 O u0 p2 c0 {7,S} {8,S}
7 O u1 p2 c0 {6,S}
8 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9 C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
"""),
    E0 = (-1081.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,180],'cm^-1')),
        HinderedRotor(inertia=(0.279325,'amu*angstrom^2'), symmetry=1, barrier=(6.42223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.445556,'amu*angstrom^2'), symmetry=1, barrier=(10.2442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.012,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3455.01,'J/mol'), sigma=(6.00748,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.66 K, Pc=36.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.41643,0.0515761,-4.23229e-05,-8.80906e-09,2.06663e-11,-130088,14.3531], Tmin=(10,'K'), Tmax=(639.073,'K')), NASAPolynomial(coeffs=[10.2735,0.0249675,-1.81505e-05,5.90961e-09,-7.1364e-13,-131298,-18.2637], Tmin=(639.073,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1081.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]OC(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1OC1F(211)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-489.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,1635.66],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2878.14,'J/mol'), sigma=(4.78691,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=449.56 K, Pc=59.54 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02952,-0.0033861,0.000100719,-1.79462e-07,1.00462e-10,-58841.1,9.38768], Tmin=(10,'K'), Tmax=(577.706,'K')), NASAPolynomial(coeffs=[2.18077,0.0257049,-1.71133e-05,5.32518e-09,-6.26769e-13,-58899.3,14.9396], Tmin=(577.706,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-489.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1OC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)[CH]F(208)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-236.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(1.00516,'amu*angstrom^2'), symmetry=1, barrier=(23.1107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.57,'J/mol'), sigma=(5.31812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.11 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19,0.034662,-3.10226e-05,1.06715e-08,-5.11771e-13,-28427.3,16.7377], Tmin=(100,'K'), Tmax=(1040.82,'K')), NASAPolynomial(coeffs=[12.0464,0.00571276,-2.17144e-06,4.35017e-10,-3.30337e-14,-30962.7,-33.5315], Tmin=(1040.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-236.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = 'F[CH]C(F)OOOC(F)(F)F(5112)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u0 p2 c0 {6,S} {7,S}
9  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
11 C u1 p0 c0 {5,S} {9,S} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1030.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,431,527,613,668,835,1199,1245,1322,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.038,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.98,'J/mol'), sigma=(5.70588,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.48 K, Pc=40.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.207389,0.0961873,-0.000118453,7.1817e-08,-1.72202e-11,-123795,35.0882], Tmin=(100,'K'), Tmax=(1015.05,'K')), NASAPolynomial(coeffs=[17.8629,0.0249765,-1.32189e-05,2.69948e-09,-1.9668e-13,-127464,-52.3599], Tmin=(1015.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1030.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFFO) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[C]C(F)(F)F(4171)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p1 c0 {4,S} {5,S}
"""),
    E0 = (-531.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,346,529,1231,1306,617,898,1187,310.428,1028.96,1029.05],'cm^-1')),
        HinderedRotor(inertia=(0.085726,'amu*angstrom^2'), symmetry=1, barrier=(5.82341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88935,0.00712406,8.594e-05,-2.10572e-07,1.46542e-10,-63912.6,11.0368], Tmin=(10,'K'), Tmax=(506.569,'K')), NASAPolynomial(coeffs=[4.91752,0.0201688,-1.53534e-05,5.20649e-09,-6.50075e-13,-64288.3,4.0956], Tmin=(506.569,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-531.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(F)OOOF(4524)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {6,S}
6  O u0 p2 c0 {3,S} {5,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u1 p0 c0 {2,S} {7,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-255.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848145,0.0699678,-8.50069e-05,4.99212e-08,-1.14684e-11,-30558.7,27.2263], Tmin=(100,'K'), Tmax=(1064.2,'K')), NASAPolynomial(coeffs=[15.4067,0.0152453,-7.87304e-06,1.59964e-09,-1.16524e-13,-33657.3,-43.9159], Tmin=(1064.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2sFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = '[CH]C(F)OOOC(F)(F)C(F)(F)F(14151)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {7,S} {8,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
12 C u0 p0 c0 {6,S} {8,S} {13,S} {14,S}
13 C u2 p0 c0 {12,S} {15,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1039.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,261,493,600,1152,1365,1422,3097,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (212.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30734,0.121973,-0.000155252,9.63404e-08,-2.35602e-11,-124811,40.1415], Tmin=(100,'K'), Tmax=(997.304,'K')), NASAPolynomial(coeffs=[21.7609,0.0294524,-1.6098e-05,3.3222e-09,-2.43227e-13,-129412,-71.0874], Tmin=(997.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1039.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCFHO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(137)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.278,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93332,-0.000263306,8.89168e-06,-1.0303e-08,3.508e-12,25853.7,4.33731], Tmin=(100,'K'), Tmax=(1056.13,'K')), NASAPolynomial(coeffs=[4.72429,0.00164127,-7.73092e-07,1.90982e-10,-1.59921e-14,25413.4,-0.815661], Tmin=(1056.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'F[CH]OOOC(F)(F)C(F)(F)F(14142)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {7,S} {8,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
12 C u1 p0 c0 {6,S} {8,S} {13,S}
13 H u0 p0 c0 {12,S}
"""),
    E0 = (-1239.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,580,1155,1237,1373,3147,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (199.029,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.637396,0.109369,-0.000153317,1.06759e-07,-2.94026e-11,-148920,36.5777], Tmin=(100,'K'), Tmax=(887.878,'K')), NASAPolynomial(coeffs=[17.4944,0.027684,-1.53199e-05,3.14552e-09,-2.28526e-13,-152140,-48.7419], Tmin=(887.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1239.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsFHHO) + radical(CsF1sHO2s)"""),
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
    label = 'F[C]C(F)OOOC(F)(F)C(F)(F)F(14152)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {13,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
14 C u2 p0 c0 {7,S} {13,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1291.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,261,493,600,1152,1365,1422,3097,90,150,250,247,323,377,431,433,1065,1274,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (230.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55032,0.129975,-0.000174097,1.13966e-07,-2.9456e-11,-155089,40.9417], Tmin=(100,'K'), Tmax=(944.263,'K')), NASAPolynomial(coeffs=[21.6519,0.0316911,-1.79735e-05,3.74335e-09,-2.74758e-13,-159471,-69.6657], Tmin=(944.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1291.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C(F)(F)C(F)(F)F(12762)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 O u1 p2 c0 {7,S}
7 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8 C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
"""),
    E0 = (-1082.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([351,323,533,609,664,892,1120,1201,193,295,551,588,656,1146,1192,1350,180],'cm^-1')),
        HinderedRotor(inertia=(0.386566,'amu*angstrom^2'), symmetry=1, barrier=(8.88792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73973,0.0181553,0.000143482,-4.22277e-07,3.34601e-10,-130149,13.0026], Tmin=(10,'K'), Tmax=(462.701,'K')), NASAPolynomial(coeffs=[7.28162,0.0284354,-2.24327e-05,7.81057e-09,-9.9484e-13,-130914,-6.08633], Tmin=(462.701,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1082.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""[O]C(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1OOC1F(210)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-446.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,598.4,1093.72,1093.74,1093.75],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3151.82,'J/mol'), sigma=(5.13508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.31 K, Pc=52.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93265,0.00400026,9.42613e-05,-1.82672e-07,1.07583e-10,-53652.7,10.67], Tmin=(10,'K'), Tmax=(549.667,'K')), NASAPolynomial(coeffs=[2.12082,0.0318962,-2.20097e-05,7.03784e-09,-8.46553e-13,-53675.7,16.3053], Tmin=(549.667,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-446.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1OOC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC=COOOC(F)(F)C(F)(F)F(14153)',
    structure = adjacencyList("""1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {7,S} {8,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
12 C u0 p0 c0 {8,S} {13,D} {14,S}
13 C u0 p0 c0 {6,S} {12,D} {15,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1283.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (212.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.894832,0.113479,-0.000141839,8.81725e-08,-2.17869e-11,-154243,39.6301], Tmin=(100,'K'), Tmax=(983.525,'K')), NASAPolynomial(coeffs=[19.2579,0.0315195,-1.6843e-05,3.44729e-09,-2.51213e-13,-158207,-57.2608], Tmin=(983.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1283.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]OOC(F)(F)C(F)(F)F(14154)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {8,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
"""),
    E0 = (-1072.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,180,180,180,625.741,880.891,880.98],'cm^-1')),
        HinderedRotor(inertia=(0.00635689,'amu*angstrom^2'), symmetry=1, barrier=(47.3665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555155,'amu*angstrom^2'), symmetry=1, barrier=(12.7641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231918,'amu*angstrom^2'), symmetry=1, barrier=(12.7629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443859,0.0811901,-0.000112916,7.5981e-08,-1.98895e-11,-128826,29.3449], Tmin=(100,'K'), Tmax=(940.238,'K')), NASAPolynomial(coeffs=[15.8983,0.0154438,-8.02989e-06,1.6133e-09,-1.16099e-13,-131732,-44.2619], Tmin=(940.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1072.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + radical(ROOJ)"""),
)

species(
    label = 'CHFCHF[Z](59)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-310.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([151,237,609,755,844,966,1147,1245,1323,1443,3181,3261],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383152,0.0294896,-2.94145e-05,1.64336e-08,-4.01759e-12,-36926.9,22.5083], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.34201,0.00821939,-3.17549e-06,5.49282e-10,-3.47434e-14,-38823.3,-13.1129], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-310.115,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCHF[Z]""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'FC=C(F)OOOC(F)(F)C(F)(F)F(14155)',
    structure = adjacencyList("""1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {13,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u0 p0 c0 {6,S} {9,S} {14,D}
14 C u0 p0 c0 {7,S} {13,D} {15,S}
15 H u0 p0 c0 {14,S}
"""),
    E0 = (-1470.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,326,540,652,719,1357,194,682,905,1196,1383,3221,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (230.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4496,0.130586,-0.000195591,1.50134e-07,-4.59895e-11,-176728,42.3439], Tmin=(100,'K'), Tmax=(798.187,'K')), NASAPolynomial(coeffs=[17.0953,0.037655,-2.09564e-05,4.28078e-09,-3.08752e-13,-179689,-42.9448], Tmin=(798.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1470.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'F[CH]C(F)OOO[C](F)C(F)(F)F(10963)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {7,S} {8,S}
10 C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {4,S} {12,S}
12 C u1 p0 c0 {5,S} {8,S} {11,S}
13 C u1 p0 c0 {6,S} {10,S} {15,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1048.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,253,522,585,610,849,1160,1215,1402,395,473,707,1436,334,575,1197,1424,3202,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (212.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67023,0.134416,-0.00019992,1.48551e-07,-4.35889e-11,-125918,42.2517], Tmin=(100,'K'), Tmax=(835.667,'K')), NASAPolynomial(coeffs=[19.4544,0.0332999,-1.84177e-05,3.75281e-09,-2.70433e-13,-129449,-55.8699], Tmin=(835.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1048.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)OOOC(F)(F)[C](F)F(14156)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {13,S}
6  F u0 p3 c0 {13,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {9,S} {11,S}
9  O u0 p2 c0 {7,S} {8,S}
10 C u0 p0 c0 {1,S} {7,S} {12,S} {14,S}
11 C u0 p0 c0 {2,S} {3,S} {8,S} {13,S}
12 C u1 p0 c0 {4,S} {10,S} {15,S}
13 C u1 p0 c0 {5,S} {6,S} {11,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-1020.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,253,525,597,667,842,1178,1324,334,575,1197,1424,3202,190,488,555,1236,1407,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (212.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90712,0.139948,-0.000211843,1.5886e-07,-4.68822e-11,-122552,42.4735], Tmin=(100,'K'), Tmax=(832.045,'K')), NASAPolynomial(coeffs=[20.4334,0.0325471,-1.822e-05,3.72063e-09,-2.68126e-13,-126270,-61.1989], Tmin=(832.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1020.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][CH]OOOC(F)(F)C(F)(F)F(14157)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {13,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {7,S} {8,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
12 C u1 p0 c0 {8,S} {13,S} {14,S}
13 C u1 p0 c0 {6,S} {12,S} {15,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1069.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,3025,407.5,1350,352.5,334,575,1197,1424,3202,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (212.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70834,0.135732,-0.000209055,1.62759e-07,-5.01174e-11,-128409,42.535], Tmin=(100,'K'), Tmax=(797.059,'K')), NASAPolynomial(coeffs=[18.408,0.0347781,-1.90678e-05,3.85158e-09,-2.75366e-13,-131615,-49.9518], Tmin=(797.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1069.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(Cs-CsOsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFHH) + radical(CCsJOO) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = '[O]OC(F)[CH]F(131)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-229.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.418102,'amu*angstrom^2'), symmetry=1, barrier=(9.61299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419586,'amu*angstrom^2'), symmetry=1, barrier=(9.64711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.65,'J/mol'), sigma=(5.51737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.71 K, Pc=45.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66565,0.0577513,-0.000101954,9.14652e-08,-3.12262e-11,-27473.1,20.2391], Tmin=(100,'K'), Tmax=(863.112,'K')), NASAPolynomial(coeffs=[7.65508,0.0169238,-8.28535e-06,1.57109e-09,-1.06542e-13,-28020.1,-4.95442], Tmin=(863.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[C](F)C(F)(F)F(11865)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 F u0 p3 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u1 p0 c0 {4,S} {5,S} {6,S}
"""),
    E0 = (-924.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.175355,'amu*angstrom^2'), symmetry=1, barrier=(9.20956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.013,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.9,'angstroms'), dipoleMoment=(1.5,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71387,0.024372,3.29209e-05,-1.25571e-07,9.2856e-11,-111156,12.4085], Tmin=(10,'K'), Tmax=(537.268,'K')), NASAPolynomial(coeffs=[6.77844,0.0217642,-1.62177e-05,5.41004e-09,-6.67368e-13,-111777,-3.18668], Tmin=(537.268,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-924.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""F[C](F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OOC(F)[CH]F(908)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-198.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180,180,860.16,860.182,860.184,860.193],'cm^-1')),
        HinderedRotor(inertia=(0.518581,'amu*angstrom^2'), symmetry=1, barrier=(11.9232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30625,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518582,'amu*angstrom^2'), symmetry=1, barrier=(11.9232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3514.65,'J/mol'), sigma=(5.71527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.98 K, Pc=42.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36451,0.0580206,-7.1795e-05,4.3694e-08,-1.03783e-11,-23756.9,24.7476], Tmin=(100,'K'), Tmax=(1033.34,'K')), NASAPolynomial(coeffs=[13.1073,0.0125657,-5.81379e-06,1.1265e-09,-8.00193e-14,-26183.8,-32.2899], Tmin=(1033.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(ROOJ) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]F(141)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-71.0739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,406,528,622,1148,1246,1368,1480,3164,3240,1663.86],'cm^-1')),
        HinderedRotor(inertia=(0.36711,'amu*angstrom^2'), symmetry=1, barrier=(8.44058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88781,0.0292781,-5.20579e-05,5.24786e-08,-2.00705e-11,-8512.84,14.1702], Tmin=(100,'K'), Tmax=(832.942,'K')), NASAPolynomial(coeffs=[3.49542,0.0155533,-7.88019e-06,1.5434e-09,-1.07577e-13,-8239.17,13.6003], Tmin=(832.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.0739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sH)"""),
)

species(
    label = 'CF3(45)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-483.19,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(68.9952,'amu')),
        NonlinearRotor(inertia=([46.5949,46.5949,90.0934],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([504.804,504.824,700.362,1092.66,1279.94,1279.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0058,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1006.05,'J/mol'), sigma=(4.32,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05587,-0.0054938,7.30062e-05,-1.34908e-07,7.857e-11,-58113,8.13986], Tmin=(10,'K'), Tmax=(570.231,'K')), NASAPolynomial(coeffs=[3.53924,0.0118582,-8.75006e-06,2.89362e-09,-3.54041e-13,-58277.3,8.38507], Tmin=(570.231,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(F)OOO[C](F)F(4537)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {8,S} {12,S}
10 C u1 p0 c0 {3,S} {4,S} {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-579.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,493,600,700,1144,1293,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0391984,0.0898682,-0.000112103,6.84585e-08,-1.64544e-11,-69516.6,34.0646], Tmin=(100,'K'), Tmax=(1015.43,'K')), NASAPolynomial(coeffs=[17.3555,0.0216562,-1.13414e-05,2.30532e-09,-1.6759e-13,-73033.4,-49.7421], Tmin=(1015.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFHO) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'F[CH][C](F)OOOC(F)(F)C(F)(F)F(14158)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {13,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u1 p0 c0 {6,S} {9,S} {14,S}
14 C u1 p0 c0 {7,S} {13,S} {15,S}
15 H u0 p0 c0 {14,S}
"""),
    E0 = (-1268.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,395,473,707,1436,334,575,1197,1424,3202,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (230.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22675,0.148013,-0.000236079,1.86821e-07,-5.80851e-11,-152348,44.8966], Tmin=(100,'K'), Tmax=(791.326,'K')), NASAPolynomial(coeffs=[20.4398,0.0334349,-1.88851e-05,3.83709e-09,-2.74464e-13,-155935,-59.1508], Tmin=(791.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1268.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = 'F[CH]C(F)OOOC(F)=C(F)F(10965)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u0 p2 c0 {6,S} {7,S}
9  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {9,S} {14,S}
11 C u0 p0 c0 {3,S} {7,S} {12,D}
12 C u0 p0 c0 {4,S} {5,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-783.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,326,540,652,719,1357,182,240,577,636,1210,1413,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (193.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.970889,0.119245,-0.000184104,1.47785e-07,-4.7446e-11,-94085,40.0395], Tmin=(100,'K'), Tmax=(762.327,'K')), NASAPolynomial(coeffs=[14.9495,0.0357006,-1.96995e-05,3.99609e-09,-2.86458e-13,-96512.1,-32.445], Tmin=(762.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-783.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[C]C(F)OOOC(F)(F)C(F)(F)F-2(14159)',
    structure = adjacencyList("""1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {13,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
14 C u0 p1 c0 {7,S} {13,S}
15 H u0 p0 c0 {13,S}
"""),
    E0 = (-1309.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,2750,2850,1437.5,1250,1305,750,350,617,898,1187,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (230.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49197,0.128537,-0.000173045,1.14374e-07,-2.98537e-11,-157347,42.055], Tmin=(100,'K'), Tmax=(935.497,'K')), NASAPolynomial(coeffs=[21.1707,0.0316393,-1.76811e-05,3.65994e-09,-2.67688e-13,-161588,-65.7687], Tmin=(935.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1309.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCFHO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'FC[C](F)OOOC(F)(F)C(F)(F)F(14145)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {14,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
14 C u1 p0 c0 {7,S} {9,S} {13,S}
15 H u0 p0 c0 {13,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1467.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (231.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7898,0.137245,-0.000194575,1.38568e-07,-3.91688e-11,-176253,42.9897], Tmin=(100,'K'), Tmax=(864.624,'K')), NASAPolynomial(coeffs=[19.9335,0.0367477,-2.02275e-05,4.13887e-09,-2.99851e-13,-180009,-58.6529], Tmin=(864.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1467.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[C](OOOC(F)C(F)F)C(F)(F)F(14160)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {13,S}
5  F u0 p3 c0 {13,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {14,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {8,S} {12,S} {15,S}
12 C u0 p0 c0 {2,S} {3,S} {11,S} {16,S}
13 C u0 p0 c0 {4,S} {5,S} {6,S} {14,S}
14 C u1 p0 c0 {7,S} {9,S} {13,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {12,S}
"""),
    E0 = (-1463.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (231.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77231,0.136297,-0.000189401,1.31422e-07,-3.61525e-11,-175818,43.1563], Tmin=(100,'K'), Tmax=(888.054,'K')), NASAPolynomial(coeffs=[20.5287,0.0358464,-1.97285e-05,4.04583e-09,-2.93839e-13,-179779,-61.7851], Tmin=(888.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1463.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'F[C](F)C(F)(F)OOOC(F)C(F)F(14161)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {13,S}
5  F u0 p3 c0 {13,S}
6  F u0 p3 c0 {14,S}
7  F u0 p3 c0 {14,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {12,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {8,S} {13,S} {15,S}
12 C u0 p0 c0 {2,S} {3,S} {9,S} {14,S}
13 C u0 p0 c0 {4,S} {5,S} {11,S} {16,S}
14 C u1 p0 c0 {6,S} {7,S} {12,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {13,S}
"""),
    E0 = (-1435.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,493,600,1152,1365,1422,3097,253,525,597,667,842,1178,1324,235,523,627,1123,1142,1372,1406,3097,190,488,555,1236,1407,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (231.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00684,0.141791,-0.00020114,1.41419e-07,-3.92795e-11,-172453,43.3703], Tmin=(100,'K'), Tmax=(881.063,'K')), NASAPolynomial(coeffs=[21.4857,0.0351344,-1.95558e-05,4.01983e-09,-2.92062e-13,-176592,-66.9924], Tmin=(881.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1435.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + radical(Csj(Cs-F1sF1sO2s)(F1s)(F1s))"""),
)

species(
    label = 'FC(F)[CH]OOOC(F)(F)C(F)(F)F(14162)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {13,S}
7  F u0 p3 c0 {13,S}
8  O u0 p2 c0 {10,S} {11,S}
9  O u0 p2 c0 {10,S} {14,S}
10 O u0 p2 c0 {8,S} {9,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {12,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
13 C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
14 C u1 p0 c0 {9,S} {13,S} {16,S}
15 H u0 p0 c0 {13,S}
16 H u0 p0 c0 {14,S}
"""),
    E0 = (-1484.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (231.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.83244,0.138387,-0.000199144,1.44696e-07,-4.17453e-11,-178366,42.8085], Tmin=(100,'K'), Tmax=(847.69,'K')), NASAPolynomial(coeffs=[19.5695,0.0373997,-2.04488e-05,4.16408e-09,-3.00552e-13,-181994,-56.9077], Tmin=(847.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1484.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + group(Cs-CsOsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CsCsFFH) + radical(CCsJOO)"""),
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
    E0 = (-736.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-356.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-111.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-304.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-362.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-417.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-710.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-503.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-709.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-597.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-313.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-285.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-334.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-649.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-460.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-656.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-481.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-400.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-394.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-130.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-438.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-427.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-656.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-602.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-588.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-623.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    products = ['[O]OC(F)(F)C(F)(F)F(13321)', 'FC1OC1F(211)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.38e+12,'s^-1','*|/',1.2), n=0, Ea=(63.953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_S;C_sec_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', 'F[CH]C(F)OOOC(F)(F)F(5112)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(215.503,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]C(F)(F)F(4171)', 'F[CH]C(F)OOOF(4524)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(12.574,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', '[CH]C(F)OOOC(F)(F)C(F)(F)F(14151)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]F(137)', 'F[CH]OOOC(F)(F)C(F)(F)F(14142)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(5)', 'F[C]C(F)OOOC(F)(F)C(F)(F)F(14152)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    products = ['[O]C(F)(F)C(F)(F)F(12762)', 'FC1OOC1F(210)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0, Ea=(89.9403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO;C_sec_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'FC=COOOC(F)(F)C(F)(F)F(14153)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(45.2717,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OOC(F)(F)C(F)(F)F(14154)', 'CHFCHF[Z](59)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.000422724,'m^3/(mol*s)'), n=2.88896, Ea=(10.4335,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.38434251561780236, var=0.3755598567717395, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0_N-4R!H->S_3OS->O_Ext-2CNS-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0_N-4R!H->S_3OS->O_Ext-2CNS-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'FC=C(F)OOOC(F)(F)C(F)(F)F(14155)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(21820,'m^3/(mol*s)'), n=0.859, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'F[CH]C(F)OOO[C](F)C(F)(F)F(10963)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -39.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[CH]C(F)OOOC(F)(F)[C](F)F(14156)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -45.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'F[CH][CH]OOOC(F)(F)C(F)(F)F(14157)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -34.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(F)(F)C(F)(F)F(12762)', '[O]OC(F)[CH]F(131)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[C](F)C(F)(F)F(11865)', '[O]OOC(F)[CH]F(908)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OC(F)(F)C(F)(F)F(13321)', '[O]C(F)[CH]F(208)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]OOC(F)(F)C(F)(F)F(14154)', 'F[CH][CH]F(141)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CF3(45)', 'F[CH]C(F)OOO[C](F)F(4537)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -16.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'F[CH][C](F)OOOC(F)(F)C(F)(F)F(14158)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R
Ea raised from -0.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['F2(78)', 'F[CH]C(F)OOOC(F)=C(F)F(10965)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CHF(40)', 'F[CH]OOOC(F)(F)C(F)(F)F(14142)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', 'F[C]C(F)OOOC(F)(F)C(F)(F)F-2(14159)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(8.03896,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    products = ['FC[C](F)OOOC(F)(F)C(F)(F)F(14145)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(34257.4,'s^-1'), n=2.45724, Ea=(144.304,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out] for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    products = ['F[C](OOOC(F)C(F)F)C(F)(F)F(14160)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(197.911,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C](F)C(F)(F)OOOC(F)C(F)F(14161)'],
    products = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(184.463,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)'],
    products = ['FC(F)[CH]OOOC(F)(F)C(F)(F)F(14162)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(177.629,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #4018',
    isomers = [
        'F[CH]C(F)OOOC(F)(F)C(F)(F)F(14144)',
    ],
    reactants = [
        ('[O]OC(F)(F)C(F)(F)F(13321)', 'FC1OC1F(211)'),
        ('[O]OC(F)(F)C(F)(F)F(13321)', '[O]C(F)[CH]F(208)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #4018',
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

