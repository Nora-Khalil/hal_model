species(
    label = '[O]C(O)O[CH]O(1091)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u1 p0 c0 {1,S} {3,S} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-355.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1682.23,1685.09,1686.89],'cm^-1')),
        HinderedRotor(inertia=(0.223055,'amu*angstrom^2'), symmetry=1, barrier=(5.12848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223722,'amu*angstrom^2'), symmetry=1, barrier=(5.1438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222205,'amu*angstrom^2'), symmetry=1, barrier=(5.10892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22278,'amu*angstrom^2'), symmetry=1, barrier=(5.12215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62902,0.070888,-0.00015662,1.71799e-07,-6.73902e-11,-42696.2,23.7109], Tmin=(100,'K'), Tmax=(859.541,'K')), NASAPolynomial(coeffs=[-1.90768,0.0406278,-2.22828e-05,4.37067e-09,-3.01166e-13,-40362.4,50.2773], Tmin=(859.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = 'O=CO(77)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-389.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,910.228,910.24],'cm^-1')),
        HinderedRotor(inertia=(0.000203447,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76624,-0.000867782,3.32891e-05,-3.95828e-08,1.39946e-11,-46830.2,7.56886], Tmin=(100,'K'), Tmax=(997.616,'K')), NASAPolynomial(coeffs=[5.52461,0.0077913,-3.35088e-06,6.86825e-10,-5.23916e-14,-47962.8,-4.82803], Tmin=(997.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]O(118)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.873,3308.85],'cm^-1')),
        HinderedRotor(inertia=(0.0103145,'amu*angstrom^2'), symmetry=1, barrier=(22.7124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75937,0.00296141,8.9037e-06,-1.3501e-08,5.39791e-12,24775.6,6.76289], Tmin=(100,'K'), Tmax=(940.438,'K')), NASAPolynomial(coeffs=[5.09115,0.00321234,-9.31657e-07,1.59608e-10,-1.15723e-14,24263.5,-0.97116], Tmin=(940.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]C([O])O(220)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {6,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (-145.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2358.11,2358.29],'cm^-1')),
        HinderedRotor(inertia=(3.03187e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79685,0.0288372,-0.000108613,1.59396e-07,-7.20139e-11,-17564.9,12.6702], Tmin=(100,'K'), Tmax=(855.716,'K')), NASAPolynomial(coeffs=[-15.1237,0.0502966,-2.88126e-05,5.74781e-09,-3.9998e-13,-11874.4,115.333], Tmin=(855.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = 'O(7)',
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
    label = 'O[CH]O[CH]O(844)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {4,S} {8,S}
4 C u1 p0 c0 {1,S} {3,S} {6,S}
5 C u1 p0 c0 {1,S} {2,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-204.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3050,390,425,1340,1360,335,370,363.296,363.301,363.304],'cm^-1')),
        HinderedRotor(inertia=(0.00127724,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114941,'amu*angstrom^2'), symmetry=1, barrier=(10.7658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114947,'amu*angstrom^2'), symmetry=1, barrier=(10.7659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114948,'amu*angstrom^2'), symmetry=1, barrier=(10.7659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18982,0.0635859,-9.99368e-05,7.54254e-08,-2.15868e-11,-24468.5,20.1187], Tmin=(100,'K'), Tmax=(943.67,'K')), NASAPolynomial(coeffs=[12.9186,0.00817872,-2.81828e-06,4.23686e-10,-2.39289e-14,-26428.7,-34.4435], Tmin=(943.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCJO)"""),
)

species(
    label = 'OC1OC(O)O1(1099)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {5,S} {9,S}
4  O u0 p2 c0 {6,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-612.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46084,0.0245792,-5.85748e-06,-2.80271e-09,9.05233e-13,-73703.9,14.1361], Tmin=(100,'K'), Tmax=(1968.47,'K')), NASAPolynomial(coeffs=[15.0224,0.0147223,-8.73769e-06,1.692e-09,-1.12559e-13,-80897.6,-56.1825], Tmin=(1968.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-612.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsOsH) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C(O)OCO(1100)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {3,S} {4,D}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-799.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5254,0.0312647,4.2097e-05,-9.72643e-08,4.49471e-11,-96036.2,21.4127], Tmin=(100,'K'), Tmax=(925.25,'K')), NASAPolynomial(coeffs=[24.0066,-0.00607479,5.6032e-06,-1.05849e-09,6.27868e-14,-102758,-99.1435], Tmin=(925.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-799.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-OsOsHH) + group(Cds-OdOsOs)"""),
)

species(
    label = 'O=COC(O)O(1095)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6  C u0 p0 c0 {1,S} {4,D} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-751.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38723,0.0304123,-1.61522e-05,3.34322e-09,-2.49347e-13,-90410.5,16.8154], Tmin=(100,'K'), Tmax=(2943.92,'K')), NASAPolynomial(coeffs=[31.3923,-0.00436415,-1.01379e-07,8.6305e-11,-4.85722e-15,-108319,-150.942], Tmin=(2943.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-751.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH)"""),
)

species(
    label = '[O][CH]O(212)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1989.6],'cm^-1')),
        HinderedRotor(inertia=(0.0936924,'amu*angstrom^2'), symmetry=1, barrier=(2.15417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40205,0.0209569,-4.96195e-05,5.95345e-08,-2.44629e-11,660.882,10.3083], Tmin=(100,'K'), Tmax=(868.742,'K')), NASAPolynomial(coeffs=[-0.343261,0.0179326,-9.40006e-06,1.81362e-09,-1.23836e-13,2076.49,32.2525], Tmin=(868.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = 'H(6)',
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
    label = '[O]C(O)OC=O(1101)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-524.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,1467.37,1467.56],'cm^-1')),
        HinderedRotor(inertia=(2.31677,'amu*angstrom^2'), symmetry=1, barrier=(53.2671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00248247,'amu*angstrom^2'), symmetry=1, barrier=(3.7923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00248446,'amu*angstrom^2'), symmetry=1, barrier=(3.79572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0428,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.10399,0.0139208,6.81949e-06,-8.61971e-09,1.73593e-12,-63140.1,11.8203], Tmin=(100,'K'), Tmax=(2142.26,'K')), NASAPolynomial(coeffs=[20.176,0.0109677,-8.75022e-06,1.71428e-09,-1.10568e-13,-75377.8,-85.8669], Tmin=(2142.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-524.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH) + radical(OCOJ)"""),
)

species(
    label = 'OH(8)',
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
    label = 'O=CO[CH]O(842)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-373.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,310.593,1117.17],'cm^-1')),
        HinderedRotor(inertia=(0.548524,'amu*angstrom^2'), symmetry=1, barrier=(38.4698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122186,'amu*angstrom^2'), symmetry=1, barrier=(8.75659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256552,'amu*angstrom^2'), symmetry=1, barrier=(22.6279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19184,0.0345799,-3.01412e-05,1.29818e-08,-2.19989e-12,-44800.1,17.8927], Tmin=(100,'K'), Tmax=(1426.66,'K')), NASAPolynomial(coeffs=[11.1079,0.00958174,-3.85827e-06,7.00192e-10,-4.77574e-14,-47344.2,-28.2905], Tmin=(1426.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-373.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-OdOsH) + radical(OCJO)"""),
)

species(
    label = 'O=C(O)O[CH]O(1102)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {2,S} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-610.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0428,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45652,0.047851,-3.92647e-05,7.56568e-09,2.73574e-12,-73321.3,21.4496], Tmin=(100,'K'), Tmax=(983.54,'K')), NASAPolynomial(coeffs=[16.0305,0.00603151,-2.10151e-06,4.16434e-10,-3.27321e-14,-77032.3,-52.9106], Tmin=(983.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-OsOsHH) + group(Cds-OdOsOs) + radical(OCJO)"""),
)

species(
    label = '[O]COC([O])O(1103)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-317.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,1307.06,1307.14,1308.5,1312.78],'cm^-1')),
        HinderedRotor(inertia=(0.153447,'amu*angstrom^2'), symmetry=1, barrier=(3.52804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156333,'amu*angstrom^2'), symmetry=1, barrier=(3.59441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154352,'amu*angstrom^2'), symmetry=1, barrier=(3.54885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25147,0.0515203,-0.00010026,1.17704e-07,-5.04503e-11,-38090.4,22.1577], Tmin=(100,'K'), Tmax=(815.791,'K')), NASAPolynomial(coeffs=[-3.05383,0.0438253,-2.41323e-05,4.84259e-09,-3.422e-13,-36103.2,53.5472], Tmin=(815.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = 'O[CH]O[C](O)O(1104)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {5,S} {8,S}
4  O u0 p2 c0 {6,S} {10,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  C u1 p0 c0 {1,S} {4,S} {7,S}
7  H u0 p0 c0 {6,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-377.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10561,0.0745098,-0.000142611,1.3441e-07,-4.72893e-11,-45333.9,25.1126], Tmin=(100,'K'), Tmax=(871.717,'K')), NASAPolynomial(coeffs=[6.9702,0.0231048,-1.20078e-05,2.29323e-09,-1.55049e-13,-45425.7,2.96252], Tmin=(871.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(Cs_P) + radical(OCJO)"""),
)

species(
    label = '[O][C](O)OCO(1105)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u1 p0 c0 {1,S} {3,S} {4,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-339.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2850,1437.5,1250,1305,750,350,360,370,350,180,180,845.482,4000],'cm^-1')),
        HinderedRotor(inertia=(0.19361,'amu*angstrom^2'), symmetry=1, barrier=(4.45148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193613,'amu*angstrom^2'), symmetry=1, barrier=(4.45154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00878286,'amu*angstrom^2'), symmetry=1, barrier=(4.45745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627609,'amu*angstrom^2'), symmetry=1, barrier=(31.8596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77708,0.0545178,-8.38589e-05,7.69477e-08,-2.88242e-11,-40730.2,24.0798], Tmin=(100,'K'), Tmax=(760.187,'K')), NASAPolynomial(coeffs=[5.67001,0.02658,-1.40242e-05,2.80576e-09,-1.99527e-13,-41106.7,7.78256], Tmin=(760.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C([O])OCO(1106)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u1 p2 c0 {6,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-317.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,180,1308.21,1308.76,1308.84,1308.94],'cm^-1')),
        HinderedRotor(inertia=(0.15454,'amu*angstrom^2'), symmetry=1, barrier=(3.55318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154687,'amu*angstrom^2'), symmetry=1, barrier=(3.55655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154614,'amu*angstrom^2'), symmetry=1, barrier=(3.55488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25147,0.0515203,-0.00010026,1.17704e-07,-5.04503e-11,-38090.4,21.4646], Tmin=(100,'K'), Tmax=(815.791,'K')), NASAPolynomial(coeffs=[-3.05383,0.0438253,-2.41323e-05,4.84259e-09,-3.422e-13,-36103.2,52.8541], Tmin=(815.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[O][CH]OC(O)O(1107)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {5,S} {10,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6  C u1 p0 c0 {1,S} {4,S} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-355.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,1682.23,1685.09,1686.89],'cm^-1')),
        HinderedRotor(inertia=(0.223055,'amu*angstrom^2'), symmetry=1, barrier=(5.12848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223722,'amu*angstrom^2'), symmetry=1, barrier=(5.1438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222205,'amu*angstrom^2'), symmetry=1, barrier=(5.10892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22278,'amu*angstrom^2'), symmetry=1, barrier=(5.12215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62902,0.070888,-0.00015662,1.71799e-07,-6.73902e-11,-42696.2,23.0178], Tmin=(100,'K'), Tmax=(859.541,'K')), NASAPolynomial(coeffs=[-1.90768,0.0406278,-2.22828e-05,4.37067e-09,-3.01166e-13,-40362.4,49.5841], Tmin=(859.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
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
    E0 = (-131.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (284.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (263.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-122.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-67.7172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-106.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-107.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-38.0722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-52.2961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-98.3903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (235.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (10.8492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (26.9191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (39.3489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (16.6905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-37.5407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(O)O[CH]O(1091)'],
    products = ['O=CO(77)', 'O=CO(77)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]O(118)', '[O]C([O])O(220)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(7)', 'O[CH]O[CH]O(844)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/O2;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(O)O[CH]O(1091)'],
    products = ['OC1OC(O)O1(1099)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(O)O[CH]O(1091)'],
    products = ['O=C(O)OCO(1100)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(O)O[CH]O(1091)'],
    products = ['O=COC(O)O(1095)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=CO(77)', '[O][CH]O(212)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(61.188,'m^3/(mol*s)'), n=1.29312, Ea=(51.8727,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.28976280499384915, var=2.1569028208455543, Tref=1000.0, N=37, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', '[O]C(O)OC=O(1101)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.07,'m^3/(mol*s)'), n=2.43, Ea=(50.0368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-4R!H-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(8)', 'O=CO[CH]O(842)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(67.9369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(6)', 'O=C(O)O[CH]O(1102)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(75.8158,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]O(212)', '[O][CH]O(212)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]COC([O])O(1103)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.41377e+07,'s^-1'), n=1.57383, Ea=(103.529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(O)O[CH]O(1091)'],
    products = ['O[CH]O[C](O)O(1104)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C](O)OCO(1105)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6e+08,'s^-1'), n=1.23, Ea=(154.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C([O])OCO(1106)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.44e+08,'s^-1'), n=1.09, Ea=(109.37,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][CH]OC(O)O(1107)'],
    products = ['[O]C(O)O[CH]O(1091)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(72901.1,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_1;O_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #582',
    isomers = [
        '[O]C(O)O[CH]O(1091)',
    ],
    reactants = [
        ('O=CO(77)', 'O=CO(77)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #582',
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

