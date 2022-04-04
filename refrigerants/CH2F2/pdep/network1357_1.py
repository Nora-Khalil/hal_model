species(
    label = '[O]C(=O)C([O])(F)C=O(3250)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {6,S}
3 O u1 p2 c0 {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-472.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,180,180,931.451,933.623,936.146,936.908],'cm^-1')),
        HinderedRotor(inertia=(2.28322,'amu*angstrom^2'), symmetry=1, barrier=(52.4957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28222,'amu*angstrom^2'), symmetry=1, barrier=(52.4728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41547,0.06574,-0.000114529,1.10893e-07,-4.17551e-11,-56702,23.9663], Tmin=(100,'K'), Tmax=(810.759,'K')), NASAPolynomial(coeffs=[4.92963,0.0306078,-1.66083e-05,3.30478e-09,-2.32421e-13,-56686.9,11.3565], Tmin=(810.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCOJ)"""),
)

species(
    label = 'CO2(14)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.56331e-09,3.13594e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88611e-07,1.21037e-10,-7.84012e-15,-14180.9,6.71043], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(=O)C([O])F(712)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-365.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([241,322,589,815,1096,1220,1302,2892,180,180,1403.11,1403.45,1403.71,1405.29],'cm^-1')),
        HinderedRotor(inertia=(0.207844,'amu*angstrom^2'), symmetry=1, barrier=(4.77874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4306.78,'J/mol'), sigma=(6.22203,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.71 K, Pc=40.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43794,0.0452692,-9.27225e-05,1.00354e-07,-3.93586e-11,-43879.1,17.7901], Tmin=(100,'K'), Tmax=(857.527,'K')), NASAPolynomial(coeffs=[0.391113,0.0277345,-1.4678e-05,2.851e-09,-1.95992e-13,-42532.3,33.1562], Tmin=(857.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + radical(O2sj(Cs-F1sCOH)) + radical(CCOJ)"""),
)

species(
    label = '[O]C(=O)C(=O)[CH]OF(3415)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,D} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {6,S} {9,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-155.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,375,552.5,462.5,1710,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728632,0.0794507,-0.000138025,1.22504e-07,-4.26717e-11,-18608.3,26.0529], Tmin=(100,'K'), Tmax=(786.409,'K')), NASAPolynomial(coeffs=[10.2988,0.0222786,-1.2772e-05,2.58672e-09,-1.8353e-13,-19850.8,-16.1479], Tmin=(786.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ) + radical(OCJC=O)"""),
)

species(
    label = '[O]C(=O)O[CH]C(=O)F(3251)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {7,D}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {2,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {3,D} {6,S}
8 C u0 p0 c0 {2,S} {4,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-527.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4301.44,'J/mol'), sigma=(6.12135,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.87 K, Pc=42.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24191,0.0591954,-6.82677e-05,3.73467e-08,-7.94166e-12,-63392,25.4372], Tmin=(100,'K'), Tmax=(1150.28,'K')), NASAPolynomial(coeffs=[14.9871,0.0113964,-5.93455e-06,1.21936e-09,-8.95679e-14,-66554,-42.7993], Tmin=(1150.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][CH]C(=O)C(=O)OF(3416)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {8,S}
6 C u0 p0 c0 {3,D} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 C u1 p0 c0 {5,S} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-173.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,375,552.5,462.5,1710,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958157,0.0738638,-0.000122779,1.09046e-07,-3.85548e-11,-20763.6,24.9253], Tmin=(100,'K'), Tmax=(773.738,'K')), NASAPolynomial(coeffs=[8.92922,0.0251908,-1.39478e-05,2.80532e-09,-1.99016e-13,-21773.6,-10.0418], Tmin=(773.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OCOJ) + radical(OCJC=O)"""),
)

species(
    label = '[O][C](F)C(=O)OC=O(3375)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,D} {7,S}
7 C u1 p0 c0 {1,S} {4,S} {6,S}
8 C u0 p0 c0 {2,S} {5,D} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-548.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48814,0.0681985,-0.00013582,1.35968e-07,-5.00421e-11,-65909.9,27.321], Tmin=(100,'K'), Tmax=(871.835,'K')), NASAPolynomial(coeffs=[3.28726,0.0289751,-1.50533e-05,2.87731e-09,-1.94937e-13,-65046.7,25.638], Tmin=(871.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdOsH) + radical(O2sj(Cs-F1sCOH)) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C([O])=C(F)C=O(2876)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {7,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {6,S} {7,D}
6 C u0 p0 c0 {4,D} {5,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-334.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([288,410,724,839,1320,2782.5,750,1395,475,1775,1000,350,440,435,1725,180.214,184.694],'cm^-1')),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1589,0.0462879,-5.61232e-05,3.81933e-08,-1.11945e-11,-40220.1,18.7663], Tmin=(100,'K'), Tmax=(805.052,'K')), NASAPolynomial(coeffs=[6.77394,0.0233574,-1.33984e-05,2.81262e-09,-2.07377e-13,-40963.2,-2.49777], Tmin=(805.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(F)([C]=O)C=O(3417)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 C u1 p0 c0 {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-272.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.10175,'amu*angstrom^2'), symmetry=1, barrier=(25.3314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08995,'amu*angstrom^2'), symmetry=1, barrier=(25.0601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14585,0.0716534,-0.000135002,1.22921e-07,-4.20125e-11,-32636.6,22.064], Tmin=(100,'K'), Tmax=(873.801,'K')), NASAPolynomial(coeffs=[8.69717,0.0177783,-9.37419e-06,1.78617e-09,-1.2023e-13,-33219.2,-9.1304], Tmin=(873.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCCJ=O)"""),
)

species(
    label = 'O=CC1(F)OOC1=O(3258)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 C u0 p0 c0 {5,D} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-609.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84824,0.0292566,3.35276e-05,-6.96678e-08,2.91393e-11,-73168.1,21.6516], Tmin=(100,'K'), Tmax=(989.872,'K')), NASAPolynomial(coeffs=[17.4218,0.010407,-4.7085e-06,1.07256e-09,-8.90212e-14,-78411,-64.2322], Tmin=(989.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-609.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(CsCCFO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cs(F)-CO) + ring(Cyclobutane)"""),
)

species(
    label = '[O]C(F)(C=O)[C]1OO1(3418)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {7,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u0 p0 c0 {5,D} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-106.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74754,0.0533679,-6.4988e-05,4.22427e-08,-1.12235e-11,-12772,25.5491], Tmin=(100,'K'), Tmax=(906.875,'K')), NASAPolynomial(coeffs=[9.23821,0.0203284,-1.03396e-05,2.0692e-09,-1.48738e-13,-14130.6,-9.85679], Tmin=(906.875,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(O2s-O2s-Cs(C-F)) + radical(C=OCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C(=O)C1(F)[CH]OO1(3419)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-226.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06351,0.0475283,-6.12528e-05,5.69601e-08,-2.37922e-11,-27189.9,21.7713], Tmin=(100,'K'), Tmax=(676.453,'K')), NASAPolynomial(coeffs=[3.80593,0.0325163,-1.75232e-05,3.57294e-09,-2.58604e-13,-27317.9,14.8426], Tmin=(676.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCOJ) + radical(CCsJOO)"""),
)

species(
    label = '[O][CH]C1(F)OOC1=O(3325)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {8,S}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 C u1 p0 c0 {5,S} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-281.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25001,0.05159,-3.95486e-05,1.18823e-08,-8.02067e-13,-33804.1,24.0148], Tmin=(100,'K'), Tmax=(1236.46,'K')), NASAPolynomial(coeffs=[16.0184,0.0144468,-7.38829e-06,1.49746e-09,-1.08644e-13,-38269,-53.6558], Tmin=(1236.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = '[O]C1(F)[CH]OOC1=O(3332)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {7,S}
3 O u0 p2 c0 {2,S} {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {6,S} {9,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-339.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81699,0.0327922,1.90581e-05,-5.39355e-08,2.39873e-11,-40745.4,22.8228], Tmin=(100,'K'), Tmax=(974.403,'K')), NASAPolynomial(coeffs=[16.1978,0.0107713,-4.02101e-06,8.38068e-10,-6.76148e-14,-45305.1,-55.1996], Tmin=(974.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(C=OCOJ) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
)

species(
    label = '[O]C(=O)C1(F)OC1[O](3420)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-421.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47536,0.031399,-1.40137e-05,1.22876e-09,1.686e-13,-50714.3,16.5879], Tmin=(100,'K'), Tmax=(2294.66,'K')), NASAPolynomial(coeffs=[24.9305,0.00624428,-5.57498e-06,1.10266e-09,-7.10317e-14,-63784.7,-111.766], Tmin=(2294.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-421.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + ring(Cs(O2)-O2s-Cs(F)) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1([O])OC1(F)C=O(3381)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8 C u0 p0 c0 {5,D} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-336.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.51985,0.0273498,-1.21442e-05,1.40689e-09,2.38929e-14,-40513.4,13.8992], Tmin=(100,'K'), Tmax=(2804.61,'K')), NASAPolynomial(coeffs=[44.8256,-0.0160082,3.48955e-06,-5.13327e-10,3.49658e-14,-68677.7,-232.024], Tmin=(2804.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsOs) + group(Cds-OdCsH) + ring(Cs(O2)-O2s-Cs(F)) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1OC(=O)C1([O])F(3371)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {7,S} {8,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-424.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91994,0.0422317,-3.55874e-05,1.63224e-08,-3.06212e-12,-51002.3,24.0102], Tmin=(100,'K'), Tmax=(1269.8,'K')), NASAPolynomial(coeffs=[9.6822,0.0177801,-6.70336e-06,1.15797e-09,-7.65595e-14,-52973.6,-15.2924], Tmin=(1269.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-424.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C=OCOJ) + radical(CCOJ) + longDistanceInteraction_cyclic(Cs(F)-CO)"""),
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
    label = '[O]C(=O)C(=O)C=O(3421)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,D}
2 O u1 p2 c0 {7,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,D} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-314.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,180,180,180,977.1,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.605031,'amu*angstrom^2'), symmetry=1, barrier=(13.9109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24951,'amu*angstrom^2'), symmetry=1, barrier=(28.7286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49266,0.0630817,-0.00011712,1.08886e-07,-3.85276e-11,-37769.1,22.7083], Tmin=(100,'K'), Tmax=(841.041,'K')), NASAPolynomial(coeffs=[7.41757,0.0191129,-1.05406e-05,2.08095e-09,-1.44404e-13,-38207.3,-1.5303], Tmin=(841.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-314.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OOJ)"""),
)

species(
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C(=O)C(=O)F(251)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
"""),
    E0 = (-483.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([286,619,818,1246,1924,377.501,378.394,379.386,379.568,380.54,2122.23],'cm^-1')),
        HinderedRotor(inertia=(0.00117727,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50218,0.0367997,-5.22213e-05,3.85636e-08,-1.15981e-11,-58127.1,16.9862], Tmin=(100,'K'), Tmax=(804.911,'K')), NASAPolynomial(coeffs=[7.32142,0.012851,-7.59242e-06,1.60056e-09,-1.1788e-13,-58902.9,-5.21811], Tmin=(804.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-483.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO) + radical(C=OC=OOJ)"""),
)

species(
    label = '[O][C]=O(145)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186753,1.03202e-05,-1.52373e-08,5.80537e-12,3804.5,8.40407], Tmin=(100,'K'), Tmax=(1021.26,'K')), NASAPolynomial(coeffs=[6.36178,0.000422737,-4.06602e-07,1.525e-10,-1.51974e-14,2816.75,-6.43921], Tmin=(1021.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = '[O]C(=O)C(=O)[C]=O(3422)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,D} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {3,D} {5,S}
7 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-201.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,1855,455,950,338.556,338.66,338.746,338.827,339.274,1595.64],'cm^-1')),
        HinderedRotor(inertia=(0.400466,'amu*angstrom^2'), symmetry=1, barrier=(32.7962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00831469,'amu*angstrom^2'), symmetry=1, barrier=(32.8129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69395,0.0563088,-9.79301e-05,8.7822e-08,-3.0846e-11,-24165,20.5411], Tmin=(100,'K'), Tmax=(796.47,'K')), NASAPolynomial(coeffs=[8.01757,0.016974,-9.5813e-06,1.92803e-09,-1.36256e-13,-24932,-7.01901], Tmin=(796.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OOJ) + radical(C=OC=OCJ=O)"""),
)

species(
    label = '[O]C(=O)C(O)(F)[C]=O(3423)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {9,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-555.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.824058,0.08151,-0.000154942,1.47653e-07,-5.31239e-11,-66750.8,23.9122], Tmin=(100,'K'), Tmax=(847.524,'K')), NASAPolynomial(coeffs=[6.98586,0.0274841,-1.51759e-05,2.98448e-09,-2.06356e-13,-66899.4,0.489744], Tmin=(847.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-555.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCOJ) + radical(C=OCCJ=O)"""),
)

species(
    label = '[O]C(F)([C]=O)C(=O)O(3424)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {7,S} {9,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-537.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689697,0.0814649,-0.000146791,1.3091e-07,-4.45281e-11,-64574.7,25.0088], Tmin=(100,'K'), Tmax=(854.003,'K')), NASAPolynomial(coeffs=[10.0678,0.0208785,-1.1111e-05,2.14804e-09,-1.46893e-13,-65569,-15.1978], Tmin=(854.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCCJ=O)"""),
)

species(
    label = '[O]C([O])=C(C=O)OF(3425)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u1 p2 c0 {8,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {7,S} {8,D}
7 C u0 p0 c0 {5,D} {6,S} {9,S}
8 C u0 p0 c0 {3,S} {4,S} {6,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-213.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,325,375,415,465,420,450,1700,1750,2782.5,750,1395,475,1775,1000,337.411,337.455,337.619],'cm^-1')),
        HinderedRotor(inertia=(0.215329,'amu*angstrom^2'), symmetry=1, barrier=(17.3902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214922,'amu*angstrom^2'), symmetry=1, barrier=(17.3883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808419,0.076884,-0.000122841,9.8297e-08,-3.11513e-11,-25520.5,23.5051], Tmin=(100,'K'), Tmax=(774.046,'K')), NASAPolynomial(coeffs=[11.9264,0.0194271,-1.14913e-05,2.38892e-09,-1.73393e-13,-27241.6,-27.2845], Tmin=(774.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=COJ)"""),
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
    E0 = (-192.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (133.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (296.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-58.3949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (285.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-67.4196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (187.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (250.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-184.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (172.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (52.8679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (49.4994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-60.0655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-90.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-56.8764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-70.6275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (44.1918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-85.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-85.5024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-181.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (96.5365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (65.7663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-117.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-65.2769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (159.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['CO2(14)', 'O=CC(=O)F(711)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(=O)C([O])F(712)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(88.9,'m^3/(mol*s)'), n=1.51, Ea=(337.662,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(=O)C(=O)[CH]OF(3415)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(173.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3251)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(134.276,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][CH]C(=O)C(=O)OF(3416)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(179.62,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O][C](F)C(=O)OC=O(3375)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.87873e+12,'s^-1'), n=0.324012, Ea=(125.251,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', '[O]C([O])=C(F)C=O(2876)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(6)', '[O]C(F)([C]=O)C=O(3417)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['O=CC1(F)OOC1=O(3258)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C(F)(C=O)[C]1OO1(3418)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(365.309,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 364.1 to 365.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C(=O)C1(F)[CH]OO1(3419)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(245.538,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 243.4 to 245.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O][CH]C1(F)OOC1=O(3325)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C1(F)[CH]OOC1=O(3332)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.24664e+09,'s^-1'), n=0.5388, Ea=(132.605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 128.9 to 132.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C(=O)C1(F)OC1[O](3420)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C1([O])OC1(F)C=O(3381)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(135.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 134.0 to 135.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C1OC(=O)C1([O])F(3371)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', '[O]C(=O)C(=O)C=O(3421)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(6.5393,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HCO(15)', '[O]C(=O)C(=O)F(251)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(86.1349,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=O(145)', 'O=CC(=O)F(711)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(86.6004,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CO2(14)', '[O][CH]C(=O)F(753)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(156.844,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C]=O(145)', '[O][CH]C(=O)F(753)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[O]C(=O)C(=O)[C]=O(3422)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(268.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C(=O)C(O)(F)[C]=O(3423)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)C([O])(F)C=O(3250)'],
    products = ['[O]C(F)([C]=O)C(=O)O(3424)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.50344e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C([O])=C(C=O)OF(3425)'],
    products = ['[O]C(=O)C([O])(F)C=O(3250)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(92.9974,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1357',
    isomers = [
        '[O]C(=O)C([O])(F)C=O(3250)',
    ],
    reactants = [
        ('CO2(14)', 'O=CC(=O)F(711)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1357',
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

