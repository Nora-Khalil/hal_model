species(
    label = '[O]OC([O])C(=O)F(3265)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {2,S}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-324.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,180,2357.11],'cm^-1')),
        HinderedRotor(inertia=(0.323391,'amu*angstrom^2'), symmetry=1, barrier=(7.4354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.323935,'amu*angstrom^2'), symmetry=1, barrier=(7.44791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73861,0.0560691,-9.73104e-05,8.85484e-08,-3.10925e-11,-38960.4,25.2293], Tmin=(100,'K'), Tmax=(838.119,'K')), NASAPolynomial(coeffs=[6.89515,0.0193857,-1.00495e-05,1.95066e-09,-1.34687e-13,-39400.7,3.79236], Tmin=(838.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-324.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(C=OCOJ) + radical(ROOJ)"""),
)

species(
    label = 'O=O(148)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,10302.3,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[O]O[CH]C(=O)F(3108)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u1 p0 c0 {2,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-215.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,180],'cm^-1')),
        HinderedRotor(inertia=(0.282257,'amu*angstrom^2'), symmetry=1, barrier=(6.48964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4772,'amu*angstrom^2'), symmetry=1, barrier=(56.9556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3770.86,'J/mol'), sigma=(5.82282,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=589.00 K, Pc=43.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06408,0.0468221,-7.59983e-05,6.61715e-08,-2.26678e-11,-25811.9,19.1781], Tmin=(100,'K'), Tmax=(821.512,'K')), NASAPolynomial(coeffs=[7.05431,0.0160915,-8.14164e-06,1.57328e-09,-1.08835e-13,-26414.8,-2.59458], Tmin=(821.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + radical(ROOJ) + radical(OCJC=O)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[O]OC([O])F(282)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-167.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,451,553,637,1069,1180,1265,1301,3056],'cm^-1')),
        HinderedRotor(inertia=(0.491221,'amu*angstrom^2'), symmetry=1, barrier=(11.2941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.73,'J/mol'), sigma=(6.011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.88 K, Pc=39.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73707,0.0289956,-3.52661e-05,2.30449e-08,-6.07484e-12,-20134.2,16.0103], Tmin=(100,'K'), Tmax=(921.449,'K')), NASAPolynomial(coeffs=[7.07196,0.0101777,-4.63275e-06,8.81523e-10,-6.15893e-14,-20933.1,-4.54836], Tmin=(921.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(ROOJ)"""),
)

species(
    label = '[O]OO[C](F)C=O(3264)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u1 p0 c0 {1,S} {2,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-153.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,376.192,378.104,378.258,758.731,1896.56],'cm^-1')),
        HinderedRotor(inertia=(0.654097,'amu*angstrom^2'), symmetry=1, barrier=(66.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662397,'amu*angstrom^2'), symmetry=1, barrier=(66.6675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119717,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3872.75,'J/mol'), sigma=(6.04376,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.91 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14412,0.0428359,-4.81855e-05,2.90158e-08,-7.13856e-12,-18348.1,22.6723], Tmin=(100,'K'), Tmax=(977.099,'K')), NASAPolynomial(coeffs=[8.55136,0.0166065,-7.91957e-06,1.54288e-09,-1.09426e-13,-19600.2,-8.09057], Tmin=(977.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C([O])C(=O)F(3353)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-304.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,180,1798.25,3800.48],'cm^-1')),
        HinderedRotor(inertia=(0.291435,'amu*angstrom^2'), symmetry=1, barrier=(6.70067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38176,0.0399151,-6.55919e-05,5.93575e-08,-2.10013e-11,-36549.3,21.4876], Tmin=(100,'K'), Tmax=(830.359,'K')), NASAPolynomial(coeffs=[5.67991,0.0160366,-8.02193e-06,1.54744e-09,-1.06926e-13,-36821.5,7.84787], Tmin=(830.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C(F)C1OOO1(4793)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {2,S} {3,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-356.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48925,0.0161545,4.79647e-05,-7.81922e-08,3.13578e-11,-42793,20.8815], Tmin=(100,'K'), Tmax=(980.296,'K')), NASAPolynomial(coeffs=[15.196,0.0068696,-2.95652e-06,7.2943e-10,-6.46605e-14,-47329.4,-50.5994], Tmin=(980.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-356.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(COCsFO) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C(F)C1OO1(3222)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-397.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57127,0.0202847,1.88708e-05,-4.79123e-08,2.25182e-11,-47757.6,17.4906], Tmin=(100,'K'), Tmax=(915.34,'K')), NASAPolynomial(coeffs=[13.4411,0.0019146,1.23703e-06,-3.00613e-10,1.83138e-14,-50967.8,-40.6548], Tmin=(915.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-397.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(COCsFO) + ring(dioxirane)"""),
)

species(
    label = 'O=C(F)C(=O)OO(4812)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {8,S}
4 O u0 p2 c0 {6,D}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {4,D} {7,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-669.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13165,0.0664156,-0.000108226,8.5239e-08,-2.58799e-11,-80382.4,20.357], Tmin=(100,'K'), Tmax=(833.061,'K')), NASAPolynomial(coeffs=[12.5592,0.0103381,-5.07955e-06,9.55615e-10,-6.46344e-14,-82244.5,-32.4358], Tmin=(833.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-669.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO)"""),
)

species(
    label = '[O]OC1OO[C]1F(4813)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u1 p2 c0 {4,S}
6 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {1,S} {3,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-20.6724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92157,0.0504869,-8.01055e-05,7.10494e-08,-2.52504e-11,-2416.08,23.0297], Tmin=(100,'K'), Tmax=(784.699,'K')), NASAPolynomial(coeffs=[6.69647,0.0201018,-1.04668e-05,2.06809e-09,-1.45571e-13,-2979.34,2.33714], Tmin=(784.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.6724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(Cs-Cs(F)-O2s-O2s) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1OOO[C]1F(4783)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {2,S} {3,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
7 C u1 p0 c0 {1,S} {3,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-55.0226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95683,0.043922,-4.42142e-05,1.87846e-08,-8.66476e-13,-6542.98,19.5129], Tmin=(100,'K'), Tmax=(801.962,'K')), NASAPolynomial(coeffs=[10.3649,0.0120354,-3.37196e-06,4.60068e-10,-2.56778e-14,-8214.76,-21.2104], Tmin=(801.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.0226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(123trioxolane) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC1OC1([O])F(4814)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u1 p2 c0 {7,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-201.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23044,0.0424349,-4.75523e-05,2.91109e-08,-7.45102e-12,-24178,20.2492], Tmin=(100,'K'), Tmax=(928.736,'K')), NASAPolynomial(coeffs=[7.669,0.0190114,-9.72091e-06,1.95473e-09,-1.41019e-13,-25188.2,-5.58665], Tmin=(928.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFOO) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(ROOJ)"""),
)

species(
    label = '[O]C1OOC1([O])F(4805)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-162.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18858,0.0455646,-7.2822e-05,7.10901e-08,-2.7891e-11,-19522.1,21.4737], Tmin=(100,'K'), Tmax=(774.545,'K')), NASAPolynomial(coeffs=[4.03762,0.025534,-1.37317e-05,2.75831e-09,-1.96253e-13,-19494.1,15.0552], Tmin=(774.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsOsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(=O)[C](O)F(4815)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u1 p0 c0 {1,S} {2,S} {7,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-407.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,280,501,1494,1531,180,180,180,180,796.596],'cm^-1')),
        HinderedRotor(inertia=(0.147438,'amu*angstrom^2'), symmetry=1, barrier=(3.38988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143018,'amu*angstrom^2'), symmetry=1, barrier=(3.28826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147367,'amu*angstrom^2'), symmetry=1, barrier=(3.38826,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36263,0.0644281,-0.000114713,9.85007e-08,-3.15155e-11,-48951.6,21.4476], Tmin=(100,'K'), Tmax=(928.042,'K')), NASAPolynomial(coeffs=[9.33421,0.0135148,-5.66416e-06,9.43267e-10,-5.70967e-14,-49718.3,-12.5745], Tmin=(928.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Os-CdOd)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + radical(C(=O)OOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'CFO(50)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC=O(704)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-118.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,524.127,524.138],'cm^-1')),
        HinderedRotor(inertia=(0.141181,'amu*angstrom^2'), symmetry=1, barrier=(27.5177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3570.08,'J/mol'), sigma=(5.61676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.64 K, Pc=45.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.11871,0.0145187,-1.8502e-07,-1.22536e-08,6.08151e-12,-14223.7,11.3177], Tmin=(100,'K'), Tmax=(990.678,'K')), NASAPolynomial(coeffs=[8.52145,0.00409734,-1.6562e-06,3.44827e-10,-2.71434e-14,-15853.2,-17.5185], Tmin=(990.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formylperoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[O]OC(=O)C(=O)F(4816)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {2,S}
6 C u0 p0 c0 {2,S} {3,D} {7,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
"""),
    E0 = (-474.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,286,619,818,1246,1924,244.171,244.447,244.461,244.464,3421.07],'cm^-1')),
        HinderedRotor(inertia=(0.265593,'amu*angstrom^2'), symmetry=1, barrier=(11.2524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265008,'amu*angstrom^2'), symmetry=1, barrier=(11.2487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36883,0.0631614,-0.000114203,9.62905e-08,-3.04021e-11,-57013.2,20.263], Tmin=(100,'K'), Tmax=(909.526,'K')), NASAPolynomial(coeffs=[10.9543,0.00870379,-4.10273e-06,7.18107e-10,-4.47651e-14,-58248.1,-22.275], Tmin=(909.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-474.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO) + radical(C(=O)OOJ)"""),
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
    label = '[O]OC(=O)[C]=O(4807)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 O u1 p2 c0 {1,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,D} {6,S}
6 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-52.7583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1855,455,950,180,180,180,3584.97],'cm^-1')),
        HinderedRotor(inertia=(0.585677,'amu*angstrom^2'), symmetry=1, barrier=(13.4659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585366,'amu*angstrom^2'), symmetry=1, barrier=(13.4587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75401,0.0533816,-9.73729e-05,8.01671e-08,-2.43577e-11,-6268.14,18.649], Tmin=(100,'K'), Tmax=(954.724,'K')), NASAPolynomial(coeffs=[10.674,0.00416036,-1.42317e-06,1.67584e-10,-5.38176e-15,-7431.35,-21.1436], Tmin=(954.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.7583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Os-CdOd)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C(=O)OOJ) + radical(OC=OCJ=O)"""),
)

species(
    label = '[O]O[C](O)C(=O)F(4817)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u1 p0 c0 {2,S} {3,S} {7,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-363.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56661,0.0593484,-0.00010261,9.08772e-08,-3.11919e-11,-43582.6,25.7545], Tmin=(100,'K'), Tmax=(828.887,'K')), NASAPolynomial(coeffs=[8.30412,0.0172705,-9.15575e-06,1.79183e-09,-1.24298e-13,-44370.9,-3.50368], Tmin=(828.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(ROOJ) + radical(Cs_P)"""),
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
    label = '[O]OC([C]=O)OF(4818)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u1 p2 c0 {2,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
7 C u1 p0 c0 {5,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-25.1227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,557,1111,1380,1390,370,380,2900,435,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.25273,'amu*angstrom^2'), symmetry=1, barrier=(5.81077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252302,'amu*angstrom^2'), symmetry=1, barrier=(5.80092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37885,'amu*angstrom^2'), symmetry=1, barrier=(31.7025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50204,0.0607918,-0.000103975,9.22933e-08,-3.20709e-11,-2937.25,25.2232], Tmin=(100,'K'), Tmax=(805.793,'K')), NASAPolynomial(coeffs=[8.3537,0.0185771,-1.01216e-05,2.01158e-09,-1.41195e-13,-3775.15,-4.70022], Tmin=(805.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.1227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C([C]=O)OOF(4819)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {5,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (27.3986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,277,555,632,1380,1390,370,380,2900,435,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.00987,'amu*angstrom^2'), symmetry=1, barrier=(23.2188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5869,0.0588786,-9.96198e-05,8.93701e-08,-3.16634e-11,3376.58,26.1499], Tmin=(100,'K'), Tmax=(788.404,'K')), NASAPolynomial(coeffs=[7.78748,0.0198436,-1.09385e-05,2.19365e-09,-1.55142e-13,2634.33,-0.796924], Tmin=(788.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.3986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCJ=O)"""),
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
    E0 = (-146.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (258.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (179.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (116.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (206.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-138.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (25.0625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-82.9676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (157.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (123.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-23.3391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (15.3818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (35.0893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-125.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-41.5182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-63.7953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-40.1843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (112.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (7.54037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (18.4259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (258.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (272.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['O=O(148)', 'O=CC(=O)F(711)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]OC([O])F(282)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(366.588,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]OO[C](F)C=O(3264)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(154.469,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', '[O]C([O])C(=O)F(3353)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', '[O]O[CH]C(=O)F(3108)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['O=C(F)C1OOO1(4793)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['O(6)', 'O=C(F)C1OO1(3222)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(171.43,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['O=C(F)C(=O)OO(4812)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['[O]OC1OO[C]1F(4813)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 303.2 to 303.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['[O]C1OOO[C]1F(4783)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(269.538,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Endocyclic
Ea raised from 268.5 to 269.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['[O]OC1OC1([O])F(4814)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(123.029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 123.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['[O]C1OOC1([O])F(4805)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(161.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 160.4 to 161.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC(=O)[C](O)F(4815)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(264.642,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)', 'O=CC(=O)F(711)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.35412e-07,'m^3/(mol*s)'), n=3.42791, Ea=(188.266,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.9589805767432776, var=2.714361599831737, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CFO(50)', '[O]OC=O(704)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(89.2092,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]OC(=O)C(=O)F(4816)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(20.9881,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[O][CH]C(=O)F(753)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(4.72618,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', '[O]OC(=O)[C]=O(4807)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(267.898,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['[O]O[C](O)C(=O)F(4817)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OC([O])C(=O)F(3265)'],
    products = ['OH(7)', '[O]C(=O)C(=O)F(251)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.83109e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC([C]=O)OF(4818)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(105.471,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C([C]=O)OOF(4819)'],
    products = ['[O]OC([O])C(=O)F(3265)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(66.6771,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1652',
    isomers = [
        '[O]OC([O])C(=O)F(3265)',
    ],
    reactants = [
        ('O=O(148)', 'O=CC(=O)F(711)'),
        ('O(6)', '[O]O[CH]C(=O)F(3108)'),
        ('O2(2)', 'O=CC(=O)F(711)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1652',
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

