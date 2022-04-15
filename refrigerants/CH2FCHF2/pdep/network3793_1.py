species(
    label = 'O=C=[C]OO[C]=C=O(14564)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {7,D}
4 O u0 p2 c0 {8,D}
5 C u1 p0 c0 {1,S} {7,D}
6 C u1 p0 c0 {2,S} {8,D}
7 C u0 p0 c0 {3,D} {5,D}
8 C u0 p0 c0 {4,D} {6,D}
"""),
    E0 = (488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1670,1700,300,440,2110,2130,495,530,650,925,180],'cm^-1')),
        HinderedRotor(inertia=(0.231777,'amu*angstrom^2'), symmetry=1, barrier=(5.32902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234396,'amu*angstrom^2'), symmetry=1, barrier=(5.38922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493917,'amu*angstrom^2'), symmetry=1, barrier=(11.3561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18536,0.0737524,-0.000148046,1.41668e-07,-5.00839e-11,58782.9,31.7147], Tmin=(100,'K'), Tmax=(876.939,'K')), NASAPolynomial(coeffs=[6.63273,0.0214421,-1.15941e-05,2.22236e-09,-1.49736e-13,58883.5,12.1705], Tmin=(876.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + missing(Cdd-CdO2d) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=C=C=O(1666)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (63.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2110,2130,495,530,650,925,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2387.04,'J/mol'), sigma=(4.99307,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=372.85 K, Pc=43.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5355,0.024033,-4.26139e-05,3.87325e-08,-1.34306e-11,7697.1,25.5769], Tmin=(100,'K'), Tmax=(857.394,'K')), NASAPolynomial(coeffs=[4.78324,0.00773532,-3.93446e-06,7.52107e-10,-5.12483e-14,7525.26,16.3243], Tmin=(857.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.731,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""OCCO(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[C]=C=O(13378)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 C u2 p0 c0 {2,D}
"""),
    E0 = (368.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,444.973],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86278,0.0119701,-1.80851e-05,1.52778e-08,-5.20063e-12,44312.6,8.89759], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.42468,0.00185394,-5.17933e-07,6.77646e-11,-3.53315e-15,43716.1,-3.69608], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(368.265,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]O[C]=C=O(8893)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
"""),
    E0 = (322.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1685,370,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.243089,'amu*angstrom^2'), symmetry=1, barrier=(5.58909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72851,0.0348827,-7.37524e-05,7.20645e-08,-2.53326e-11,38801.1,18.3764], Tmin=(100,'K'), Tmax=(910.539,'K')), NASAPolynomial(coeffs=[4.31527,0.0105752,-5.14874e-06,9.24356e-10,-5.88728e-14,39230.8,14.8163], Tmin=(910.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C=C1OOC1=C=O(14565)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {7,D}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {6,S} {7,D}
6 C u0 p0 c0 {2,S} {5,S} {8,D}
7 C u0 p0 c0 {3,D} {5,D}
8 C u0 p0 c0 {4,D} {6,D}
"""),
    E0 = (117.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51832,0.0457764,-3.35843e-05,3.16874e-09,3.65095e-12,14255.5,21.4967], Tmin=(100,'K'), Tmax=(1013.61,'K')), NASAPolynomial(coeffs=[15.8619,0.0071783,-3.11092e-06,6.51609e-10,-5.08025e-14,10422.8,-52.4596], Tmin=(1013.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)(Cds-Cdd-O2d)O2s) + group(Cds-(Cdd-O2d)(Cds-Cdd-O2d)O2s) + missing(Cdd-CdO2d) + missing(Cdd-CdO2d) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C=[C]OOC1=[C]O1(14567)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u1 p0 c0 {1,S} {5,D}
7 C u1 p0 c0 {3,S} {8,D}
8 C u0 p0 c0 {4,D} {7,D}
"""),
    E0 = (680.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12015,0.0770572,-0.000159688,1.56524e-07,-5.63484e-11,81940.5,30.8231], Tmin=(100,'K'), Tmax=(872.792,'K')), NASAPolynomial(coeffs=[5.67326,0.0240603,-1.33873e-05,2.59605e-09,-1.76182e-13,82369.5,16.4872], Tmin=(872.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + ring(Cyclopropene) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=C=C1O[C]=[C]OO1(14568)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {2,S} {8,D}
6 C u1 p0 c0 {1,S} {7,D}
7 C u1 p0 c0 {3,S} {6,D}
8 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (420.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85566,0.0260667,4.3881e-05,-8.68141e-08,3.7182e-11,50717.8,24.1595], Tmin=(100,'K'), Tmax=(964.21,'K')), NASAPolynomial(coeffs=[20.3228,0.00196844,-3.21817e-07,2.31047e-10,-3.18549e-14,44715.5,-76.919], Tmin=(964.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + ring(Cyclohexane) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O]C#C[O](9818)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
"""),
    E0 = (9.55717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,344.097,344.099,344.1],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09811,0.0228497,-4.17927e-05,3.86217e-08,-1.34608e-11,1179.06,7.869], Tmin=(100,'K'), Tmax=(871.084,'K')), NASAPolynomial(coeffs=[4.9831,0.00749284,-3.80902e-06,7.20359e-10,-4.85535e-14,1104.9,0.494495], Tmin=(871.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.55717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""OCCO""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[C]#COO[C]=C=O(14569)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 C u0 p0 c0 {2,S} {7,T}
7 C u1 p0 c0 {6,T}
"""),
    E0 = (844.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,2120,512.5,787.5,2175,525,180],'cm^-1')),
        HinderedRotor(inertia=(0.551846,'amu*angstrom^2'), symmetry=1, barrier=(12.688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.46079,'amu*angstrom^2'), symmetry=1, barrier=(56.5784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.45981,'amu*angstrom^2'), symmetry=1, barrier=(56.5559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31409,0.0695071,-0.000138723,1.2978e-07,-4.47767e-11,101625,25.812], Tmin=(100,'K'), Tmax=(890.391,'K')), NASAPolynomial(coeffs=[7.50944,0.0171043,-9.04875e-06,1.6962e-09,-1.11947e-13,101495,2.11168], Tmin=(890.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(844.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(Ct-CtOs) + missing(Cdd-CdO2d) + group(Ct-CtH) + radical(C=CJO) + radical(Acetyl)"""),
)

species(
    label = '[O]C1=[C]OOC1=C=O(14570)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {6,S} {8,D}
6 C u0 p0 c0 {3,S} {5,S} {7,D}
7 C u1 p0 c0 {2,S} {6,D}
8 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (238.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81069,0.0399676,-2.29998e-05,-8.56338e-09,9.01297e-12,28769,24.8856], Tmin=(100,'K'), Tmax=(925.223,'K')), NASAPolynomial(coeffs=[14.815,0.00494357,-5.83537e-07,4.67251e-11,-4.36926e-15,25455.3,-41.7451], Tmin=(925.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(C=CJO)"""),
)

species(
    label = '[C]1=[C]OOC#COO1(14571)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u0 p2 c0 {4,S} {8,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u1 p0 c0 {1,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 C u0 p0 c0 {2,S} {8,T}
8 C u0 p0 c0 {3,S} {7,T}
"""),
    E0 = (808.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00826,0.0725882,-0.000123189,1.09243e-07,-3.86488e-11,97293.4,15.1725], Tmin=(100,'K'), Tmax=(743.689,'K')), NASAPolynomial(coeffs=[9.54735,0.0223245,-1.30636e-05,2.68487e-09,-1.92842e-13,96143.2,-22.6889], Tmin=(743.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtOs) + ring(Ring) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]C1=C([C]=O)OO1(14572)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {7,D}
4 O u0 p2 c0 {8,D}
5 C u0 p0 c0 {1,S} {6,D} {7,S}
6 C u0 p0 c0 {2,S} {5,D} {8,S}
7 C u1 p0 c0 {3,D} {5,S}
8 C u1 p0 c0 {4,D} {6,S}
"""),
    E0 = (222.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646809,0.0820618,-0.000145082,1.25217e-07,-4.26878e-11,26891.6,22.4571], Tmin=(100,'K'), Tmax=(720.741,'K')), NASAPolynomial(coeffs=[12.1058,0.0184589,-1.26969e-05,2.74935e-09,-2.02945e-13,25240,-29.0718], Tmin=(720.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(C=CCJ=O) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=C1[C]OOC#CO1(14573)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {8,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {1,S} {4,D} {6,S}
6 C u2 p0 c0 {2,S} {5,S}
7 C u0 p0 c0 {1,S} {8,T}
8 C u0 p0 c0 {3,S} {7,T}
"""),
    E0 = (489.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46404,0.0478788,-3.41755e-05,7.69461e-09,-9.56137e-14,58956,15.6051], Tmin=(100,'K'), Tmax=(1449.97,'K')), NASAPolynomial(coeffs=[19.7896,0.0100996,-8.30863e-06,1.87792e-09,-1.40391e-13,52298.8,-84.2447], Tmin=(1449.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCt) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Ct-CtOs) + group(Ct-CtOs) + ring(cycloheptyne) + radical(CH2_triplet)"""),
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
    E0 = (157.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (360.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (165.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (353.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (223.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (157.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (157.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (756.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (216.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (478.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (237.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (288.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['O=C=C=O(1666)', 'O=C=C=O(1666)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=C=O(13378)', '[O]O[C]=C=O(8893)'],
    products = ['O=C=[C]OO[C]=C=O(14564)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['O=C=C1OOC1=C=O(14565)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['O=C=[C]OOC1=[C]O1(14567)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.262e+12,'s^-1'), n=0.216, Ea=(196.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3;multiplebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['O=C=C1O[C]=[C]OO1(14568)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.65022e+10,'s^-1'), n=0.274726, Ea=(66.2661,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;multiplebond_intra;radadd_intra]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C=C=O(1666)', '[O]C#C[O](9818)'],
    products = ['O=C=[C]OO[C]=C=O(14564)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.64086e-07,'m^3/(mol*s)'), n=3.71185, Ea=(414.712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0
Ea raised from 410.1 to 414.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C#C[O](9818)', '[O]C#C[O](9818)'],
    products = ['O=C=[C]OO[C]=C=O(14564)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(468.886,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0
Ea raised from 463.3 to 468.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(6)', '[C]#COO[C]=C=O(14569)'],
    products = ['O=C=[C]OO[C]=C=O(14564)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Ct_rad;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['[O]C1=[C]OOC1=C=O(14570)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.94e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_T;triplebond_intra;radadd_intra]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['[C]1=[C]OOC#COO1(14571)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.94988e+10,'s^-1'), n=0.402453, Ea=(320.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;multiplebond_intra;radadd_intra] for rate rule [R8;multiplebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['O=[C]C1=C([C]=O)OO1(14572)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_T;triplebond_intra;radadd_intra]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C=[C]OO[C]=C=O(14564)'],
    products = ['O=C1[C]OOC#CO1(14573)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.14e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra;radadd_intra] for rate rule [R7_linear;triplebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #3793',
    isomers = [
        'O=C=[C]OO[C]=C=O(14564)',
    ],
    reactants = [
        ('O=C=C=O(1666)', 'O=C=C=O(1666)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3793',
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

