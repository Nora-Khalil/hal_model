species(
    label = '[CH2]C(F)(F)O[CH]C(=O)F(4327)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  C u1 p0 c0 {6,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-657.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,611,648,830,1210,1753,272.892,273.358,273.81,273.886],'cm^-1')),
        HinderedRotor(inertia=(0.00225198,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293682,'amu*angstrom^2'), symmetry=1, barrier=(15.7149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152524,'amu*angstrom^2'), symmetry=1, barrier=(8.05506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.52836,'amu*angstrom^2'), symmetry=1, barrier=(28.18,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162052,0.0896937,-0.000125478,8.84016e-08,-2.46017e-11,-78985,29.1774], Tmin=(100,'K'), Tmax=(880.27,'K')), NASAPolynomial(coeffs=[14.8311,0.0230362,-1.18911e-05,2.37686e-09,-1.70125e-13,-81567.5,-39.7214], Tmin=(880.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-657.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(COCsFO) + radical(CCsJOCs) + radical(Csj(Cs-F1sF1sO2s)(H)(H))"""),
)

species(
    label = 'O=CC(=O)F(4234)',
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
    label = 'CH2CF2(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.616,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C(=O)F(537)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-5.0725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([486,617,768,1157,1926,180,1655.05,1655.52],'cm^-1')),
        HinderedRotor(inertia=(0.0191536,'amu*angstrom^2'), symmetry=1, barrier=(5.31745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3044.49,'J/mol'), sigma=(4.92747,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.54 K, Pc=57.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32766,0.0152939,-1.10758e-05,3.91582e-09,-5.6142e-13,-586.537,12.101], Tmin=(100,'K'), Tmax=(1580.41,'K')), NASAPolynomial(coeffs=[6.53552,0.00717479,-3.36979e-06,6.65147e-10,-4.72043e-14,-1600.48,-4.84321], Tmin=(1580.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])(F)F(162)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-268.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.181121,'amu*angstrom^2'), symmetry=1, barrier=(4.16432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3215.71,'J/mol'), sigma=(5.6231,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.29 K, Pc=41.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24484,0.0356831,-3.90797e-05,2.09126e-08,-4.3132e-12,-32216,16.3892], Tmin=(100,'K'), Tmax=(1195.23,'K')), NASAPolynomial(coeffs=[10.9835,0.00643772,-2.37702e-06,4.40782e-10,-3.12139e-14,-34304.9,-27.3284], Tmin=(1195.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(H)(H))"""),
)

species(
    label = 'CH2(T)(18)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(F)[CH]O[C](F)F(4080)',
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
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,493,600,700,1144,1293,244.153,244.166,244.198,244.198],'cm^-1')),
        HinderedRotor(inertia=(0.237006,'amu*angstrom^2'), symmetry=1, barrier=(10.0272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00235805,'amu*angstrom^2'), symmetry=1, barrier=(10.0269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645282,'amu*angstrom^2'), symmetry=1, barrier=(27.298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888445,0.0741621,-0.000124068,1.048e-07,-3.46432e-11,-73370.2,26.6211], Tmin=(100,'K'), Tmax=(800.21,'K')), NASAPolynomial(coeffs=[11.1627,0.0179605,-9.63843e-06,1.90235e-09,-1.32985e-13,-74859.4,-19.6874], Tmin=(800.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFHO) + group(COCsFO) + radical(CCsJOCs) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O=C(F)C1CC(F)(F)O1(4335)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
9  C u0 p0 c0 {3,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-949.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34606,0.0459707,-8.08012e-06,-2.35836e-08,1.19603e-11,-114076,27.0191], Tmin=(100,'K'), Tmax=(1040.66,'K')), NASAPolynomial(coeffs=[15.3665,0.0191755,-8.51305e-06,1.71335e-09,-1.27355e-13,-118461,-48.2295], Tmin=(1040.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-949.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(COCsFO) + ring(Cs-Cs(F)(F)-O2s-Cs)"""),
)

species(
    label = '[CH2]C(F)(F)OC1O[C]1F(4460)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
8  C u1 p0 c0 {3,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-507.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895328,0.0719558,-7.73438e-05,4.25268e-08,-9.51266e-12,-60958.6,26.035], Tmin=(100,'K'), Tmax=(1067.47,'K')), NASAPolynomial(coeffs=[12.987,0.0266461,-1.3675e-05,2.7637e-09,-2.00193e-13,-63540.1,-33.0898], Tmin=(1067.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-507.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsCFFO) + group(Cs-CsHHH) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-F1sF1sO2s)(H)(H))"""),
)

species(
    label = 'F[C]1[CH]OC(F)(F)CO1(4378)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {4,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-650.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.102671,0.0814854,-0.000100001,6.24823e-08,-1.43806e-11,-78033.4,22.885], Tmin=(100,'K'), Tmax=(1284.83,'K')), NASAPolynomial(coeffs=[15.4966,0.0144445,-1.61401e-07,-5.14412e-10,5.50461e-14,-80516.8,-50.3473], Tmin=(1284.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFFO) + group(CsCFHO) + group(Cs-CsOsHH) + ring(1,4-Dioxane) + radical(CsCsF1sO2s) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1(F)[CH]OC(F)(F)C1(4417)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
9  C u1 p0 c0 {4,S} {7,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-660.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0402667,0.0729988,-7.6183e-05,3.96309e-08,-7.77711e-12,-79332,23.4465], Tmin=(100,'K'), Tmax=(1422.69,'K')), NASAPolynomial(coeffs=[18.29,0.0121164,-1.90027e-06,9.32339e-11,1.60141e-15,-83556.1,-67.6277], Tmin=(1422.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-660.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCFFO) + ring(Tetrahydrofuran) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2][C](F)F(141)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25818,0.0186259,-1.56413e-05,8.22305e-09,-2.04579e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.66,'K')), NASAPolynomial(coeffs=[4.47009,0.0131647,-6.4128e-06,1.29206e-09,-9.37476e-14,-13305.9,7.85289], Tmin=(887.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
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
    label = 'C=C(F)O[CH]C(=O)F(4461)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u1 p0 c0 {3,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {3,S} {8,D}
7  C u0 p0 c0 {2,S} {4,D} {5,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-497.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,326,540,652,719,1357,611,648,830,1210,1753,2950,3100,1380,975,1025,1650,271.899,271.899,271.899,271.899],'cm^-1')),
        HinderedRotor(inertia=(0.356002,'amu*angstrom^2'), symmetry=1, barrier=(18.6765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356002,'amu*angstrom^2'), symmetry=1, barrier=(18.6765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356002,'amu*angstrom^2'), symmetry=1, barrier=(18.6765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448829,0.0744385,-8.73475e-05,4.878e-08,-1.0506e-11,-59688.4,25.4496], Tmin=(100,'K'), Tmax=(1143.24,'K')), NASAPolynomial(coeffs=[18.2592,0.0121225,-5.58447e-06,1.10046e-09,-7.95275e-14,-63760.6,-62.8592], Tmin=(1143.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-497.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(CdCFO) + group(COCsFO) + group(Cds-CdsHH) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][CH]C(=O)F(398)',
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
    label = '[CH2]C(F)(F)OC=C=O(4462)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-457.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04107,'amu*angstrom^2'), symmetry=1, barrier=(23.9362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0367,'amu*angstrom^2'), symmetry=1, barrier=(23.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04018,'amu*angstrom^2'), symmetry=1, barrier=(23.9158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0661167,0.0782139,-9.16824e-05,4.98439e-08,-1.02068e-11,-54827,25.7639], Tmin=(100,'K'), Tmax=(1270.85,'K')), NASAPolynomial(coeffs=[22.2141,0.0046812,-8.70749e-07,9.68119e-11,-5.7827e-15,-60214.9,-85.9832], Tmin=(1270.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-457.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sO2s)(H)(H))"""),
)

species(
    label = '[CH2][C](F)OC(F)C(=O)F(4463)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
7  C u1 p0 c0 {2,S} {4,S} {9,S}
8  C u0 p0 c0 {3,S} {5,D} {6,S}
9  C u1 p0 c0 {7,S} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-641.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,395,473,707,1436,486,617,768,1157,1926,3000,3100,440,815,1455,1000,226.403,226.426,226.468,1692.7],'cm^-1')),
        HinderedRotor(inertia=(0.240786,'amu*angstrom^2'), symmetry=1, barrier=(8.76388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240837,'amu*angstrom^2'), symmetry=1, barrier=(8.76409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939759,'amu*angstrom^2'), symmetry=1, barrier=(34.1941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940089,'amu*angstrom^2'), symmetry=1, barrier=(34.1946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.134001,0.0936217,-0.000154005,1.35541e-07,-4.73453e-11,-76974.2,31.0991], Tmin=(100,'K'), Tmax=(786.181,'K')), NASAPolynomial(coeffs=[10.3656,0.0316884,-1.69954e-05,3.38104e-09,-2.38362e-13,-78277.7,-13.86], Tmin=(786.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-641.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(COCsFO) + radical(Csj(Cs-HHH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sH)(H)(H))"""),
)

species(
    label = 'O=C(F)[CH]O[C](F)CF(4464)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {4,S} {6,S}
8  C u1 p0 c0 {4,S} {9,S} {12,S}
9  C u0 p0 c0 {3,S} {5,D} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-618.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,3025,407.5,1350,352.5,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.142031,0.0972686,-0.000150969,1.18304e-07,-3.64043e-11,-74243.6,29.894], Tmin=(100,'K'), Tmax=(800.282,'K')), NASAPolynomial(coeffs=[14.7218,0.0229721,-1.17056e-05,2.28728e-09,-1.60059e-13,-76622.6,-38.5033], Tmin=(800.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-618.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-(Cds-O2d)OsHH) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO) + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs)) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(F)(F)OC(F)[C]=O(4449)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
8  C u1 p0 c0 {6,S} {11,S} {12,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-654.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,355,410,600,1181,1341,1420,3056,3000,3100,440,815,1455,1000,1855,455,950,229.859,245.579,251.537],'cm^-1')),
        HinderedRotor(inertia=(0.208336,'amu*angstrom^2'), symmetry=1, barrier=(8.60448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201643,'amu*angstrom^2'), symmetry=1, barrier=(8.58468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203437,'amu*angstrom^2'), symmetry=1, barrier=(8.59272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.898623,'amu*angstrom^2'), symmetry=1, barrier=(37.1591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0739084,0.0978523,-0.000162878,1.39964e-07,-4.7065e-11,-78567.4,29.4371], Tmin=(100,'K'), Tmax=(815.99,'K')), NASAPolynomial(coeffs=[12.1093,0.0275951,-1.43614e-05,2.80302e-09,-1.94726e-13,-80205,-24.7133], Tmin=(815.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-654.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sO2s)(H)(H)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O][C]=COC(F)(F)CF(4465)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-642.453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,528,1116,1182,1331,1402,1494,3075,3110,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23757,'amu*angstrom^2'), symmetry=1, barrier=(28.4542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23862,'amu*angstrom^2'), symmetry=1, barrier=(28.4783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23857,'amu*angstrom^2'), symmetry=1, barrier=(28.4773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.947922,0.0886825,-9.97571e-05,5.08996e-08,-9.61863e-12,-77073.9,32.6349], Tmin=(100,'K'), Tmax=(1463.97,'K')), NASAPolynomial(coeffs=[26.8751,0.00141959,1.17247e-06,-3.07549e-10,2.17205e-14,-84015.5,-108.085], Tmin=(1463.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-642.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(F)(F)C([O])(F)C=O(4326)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
8  C u1 p0 c0 {7,S} {11,S} {12,S}
9  C u0 p0 c0 {5,D} {6,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-598.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.900228,'amu*angstrom^2'), symmetry=1, barrier=(20.698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890917,'amu*angstrom^2'), symmetry=1, barrier=(20.4839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08069,'amu*angstrom^2'), symmetry=1, barrier=(24.8473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.177957,0.092345,-0.000138056,1.02839e-07,-2.82974e-11,-71877.8,28.806], Tmin=(100,'K'), Tmax=(650.903,'K')), NASAPolynomial(coeffs=[12.6366,0.0271741,-1.41213e-05,2.78907e-09,-1.96585e-13,-73741,-27.8037], Tmin=(650.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C(F)(F)OC=[C]F(4466)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {8,D} {11,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-335.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.915153,'amu*angstrom^2'), symmetry=1, barrier=(21.0412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916416,'amu*angstrom^2'), symmetry=1, barrier=(21.0702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916482,'amu*angstrom^2'), symmetry=1, barrier=(21.0717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.409344,0.0788416,-9.92064e-05,6.03824e-08,-1.42682e-11,-40266.6,25.2225], Tmin=(100,'K'), Tmax=(1040.7,'K')), NASAPolynomial(coeffs=[17.0269,0.0149704,-7.14552e-06,1.40817e-09,-1.01083e-13,-43725.3,-55.6103], Tmin=(1040.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sF1sO2s)(H)(H)) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'FC1=COC(F)(F)CO1(4387)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-890.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.158981,0.0688453,-6.59663e-05,3.11002e-08,-5.45087e-12,-106985,19.7569], Tmin=(100,'K'), Tmax=(1657.4,'K')), NASAPolynomial(coeffs=[18.0781,0.0109694,-1.04096e-06,-6.14411e-11,1.0692e-14,-111126,-71.697], Tmin=(1657.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFFO) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(23dihydro14dioxin)"""),
)

species(
    label = '[O][C](F)C1CC(F)(F)O1(4467)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
9  C u1 p0 c0 {3,S} {5,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-579.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26992,0.0640491,-6.76738e-05,3.89786e-08,-9.34383e-12,-69611.1,26.469], Tmin=(100,'K'), Tmax=(990.76,'K')), NASAPolynomial(coeffs=[10.0853,0.028458,-1.37884e-05,2.71938e-09,-1.94351e-13,-71357.9,-15.9783], Tmin=(990.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCFHO) + ring(Cs-Cs(F)(F)-O2s-Cs) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    label = '[CH2]C(F)(F)O[C]=C=O(4455)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {8,D}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-217.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,263.059,263.072,263.11],'cm^-1')),
        HinderedRotor(inertia=(0.366372,'amu*angstrom^2'), symmetry=1, barrier=(17.9958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366442,'amu*angstrom^2'), symmetry=1, barrier=(17.9954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366314,'amu*angstrom^2'), symmetry=1, barrier=(17.9957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613233,0.0744049,-9.89851e-05,6.31285e-08,-1.55152e-11,-26026.3,26.4574], Tmin=(100,'K'), Tmax=(1005.36,'K')), NASAPolynomial(coeffs=[16.3894,0.0116368,-5.33495e-06,1.02801e-09,-7.28209e-14,-29198.4,-49.7377], Tmin=(1005.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sO2s)(H)(H)) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(F)(F)O[C]=C(O)F(4468)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u1 p0 c0 {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-557.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,293,496,537,1218,1685,370,180,180.02,180.046,774.77],'cm^-1')),
        HinderedRotor(inertia=(0.0533782,'amu*angstrom^2'), symmetry=1, barrier=(22.8035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.590118,'amu*angstrom^2'), symmetry=1, barrier=(13.574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594538,'amu*angstrom^2'), symmetry=1, barrier=(22.8044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.590516,'amu*angstrom^2'), symmetry=1, barrier=(13.5797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0307002,0.0906289,-0.000119642,7.66243e-08,-1.90976e-11,-66883.9,29.9863], Tmin=(100,'K'), Tmax=(986.445,'K')), NASAPolynomial(coeffs=[17.9316,0.0177925,-8.88684e-06,1.77285e-09,-1.27615e-13,-70427.7,-56.4263], Tmin=(986.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-557.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(Cs-F1sF1sO2s)(H)(H)) + radical(C=CJO)"""),
)

species(
    label = 'CC(F)(F)O[C]=C([O])F(4469)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-631.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,2750,2800,2850,1350,1500,750,1050,1375,1000,293,496,537,1218,1685,370,291.931,291.932,291.932,291.933,604.705],'cm^-1')),
        HinderedRotor(inertia=(0.0104484,'amu*angstrom^2'), symmetry=1, barrier=(20.2046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334087,'amu*angstrom^2'), symmetry=1, barrier=(20.2046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.33409,'amu*angstrom^2'), symmetry=1, barrier=(20.2046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.585202,0.079018,-9.34928e-05,5.58241e-08,-1.33567e-11,-75882.7,27.7472], Tmin=(100,'K'), Tmax=(1010.38,'K')), NASAPolynomial(coeffs=[14.3095,0.0246864,-1.28351e-05,2.6063e-09,-1.89289e-13,-78656.1,-38.6067], Tmin=(1010.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-631.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C(F)O[CH][C](F)OF(4470)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {3,S} {8,S}
6  C u1 p0 c0 {4,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {4,S} {9,D}
8  C u1 p0 c0 {1,S} {5,S} {6,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-233.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,326,540,652,719,1357,395,473,707,1436,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.550857,0.106801,-0.000169474,1.31985e-07,-4.00536e-11,-27897,29.7589], Tmin=(100,'K'), Tmax=(813.151,'K')), NASAPolynomial(coeffs=[16.8873,0.0210143,-1.12147e-05,2.22623e-09,-1.56983e-13,-30732.8,-50.7623], Tmin=(813.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-CsOsHH) + group(CsCFHO) + group(CdCFO) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH2]C(F)(F)OC=[C]OF(4471)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {3,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u1 p0 c0 {6,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {9,D} {12,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-241.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,283.025,283.032,283.034,283.036],'cm^-1')),
        HinderedRotor(inertia=(0.357886,'amu*angstrom^2'), symmetry=1, barrier=(20.3437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357892,'amu*angstrom^2'), symmetry=1, barrier=(20.3437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357843,'amu*angstrom^2'), symmetry=1, barrier=(20.3437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357871,'amu*angstrom^2'), symmetry=1, barrier=(20.3438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.387646,0.0915342,-0.000115377,6.84278e-08,-1.54827e-11,-28832.9,31.7082], Tmin=(100,'K'), Tmax=(1096.49,'K')), NASAPolynomial(coeffs=[21.9168,0.0101684,-4.06938e-06,7.53703e-10,-5.31349e-14,-33724.2,-77.9528], Tmin=(1096.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFFO) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Csj(Cs-F1sF1sO2s)(H)(H)) + radical(C=CJO)"""),
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
    E0 = (-275.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (108.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (152.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-267.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-91.9002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-187.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-193.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-197.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (7.74122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-193.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (38.1495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (58.9488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-39.2352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-101.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-85.5843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-68.7329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-83.4164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (289.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-267.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-197.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (174.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-24.4788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-156.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (198.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (193.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['O=CC(=O)F(4234)', 'CH2CF2(57)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=O)F(537)', '[CH2]C([O])(F)F(162)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(18)', 'O=C(F)[CH]O[C](F)F(4080)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['O=C(F)C1CC(F)(F)O1(4335)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['[CH2]C(F)(F)OC1O[C]1F(4460)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['F[C]1[CH]OC(F)(F)CO1(4378)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['[O]C1(F)[CH]OC(F)(F)C1(4417)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.93454e+07,'s^-1'), n=0.812267, Ea=(82.1501,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=CC(=O)F(4234)', '[CH2][C](F)F(141)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(61.188,'m^3/(mol*s)'), n=1.29312, Ea=(11.9067,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.28976280499384915, var=2.1569028208455543, Tref=1000.0, N=37, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'C=C(F)O[CH]C(=O)F(4461)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(49.7427,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][CH]C(=O)F(398)', 'CH2CF2(57)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[CH2]C(F)(F)OC=C=O(4462)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.9385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][CH]C(=O)F(398)', '[CH2][C](F)F(141)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C](F)OC(F)C(=O)F(4463)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(219.379,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C(F)[CH]O[C](F)CF(4464)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(134.386,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(F)(F)OC(F)[C]=O(4449)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(186.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=COC(F)(F)CF(4465)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(191.246,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(F)(F)C([O])(F)C=O(4326)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(132.816,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(6)', '[CH2]C(F)(F)OC=[C]F(4466)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['FC1=COC(F)(F)CO1(4387)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Ypri_rad_out] + [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    products = ['[O][C](F)C1CC(F)(F)O1(4467)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(78.26,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 74.7 to 78.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[CH2]C(F)(F)O[C]=C=O(4455)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(290.474,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(F)(F)O[C]=C(O)F(4468)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC(F)(F)O[C]=C([O])F(4469)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.74832e+07,'s^-1'), n=1.435, Ea=(93,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS_OCs;Y_rad_out;Cs_H_out_2H] + [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(F)O[CH][C](F)OF(4470)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(49.5471,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(F)(F)OC=[C]OF(4471)'],
    products = ['[CH2]C(F)(F)O[CH]C(=O)F(4327)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(52.0458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1449',
    isomers = [
        '[CH2]C(F)(F)O[CH]C(=O)F(4327)',
    ],
    reactants = [
        ('O=CC(=O)F(4234)', 'CH2CF2(57)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1449',
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

