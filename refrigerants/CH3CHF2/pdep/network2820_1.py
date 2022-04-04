species(
    label = 'O=C(F)OC(O)O(7264)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,D}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-974.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,482,664,788,1296,1923,180,832.438,837.788],'cm^-1')),
        HinderedRotor(inertia=(0.164646,'amu*angstrom^2'), symmetry=1, barrier=(3.78555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167392,'amu*angstrom^2'), symmetry=1, barrier=(3.84868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169021,'amu*angstrom^2'), symmetry=1, barrier=(3.88612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16384,'amu*angstrom^2'), symmetry=1, barrier=(3.767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68105,0.0568874,-8.83044e-05,8.11301e-08,-3.04016e-11,-117110,23.7156], Tmin=(100,'K'), Tmax=(757.681,'K')), NASAPolynomial(coeffs=[5.84978,0.0272479,-1.4518e-05,2.91352e-09,-2.07497e-13,-117522,6.20626], Tmin=(757.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-974.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(COFOO)"""),
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
    label = 'O=CO(1075)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76621,-0.000867468,3.3288e-05,-3.95813e-08,1.3994e-11,-46830.2,7.56895], Tmin=(100,'K'), Tmax=(997.627,'K')), NASAPolynomial(coeffs=[5.52472,0.00779112,-3.35078e-06,6.868e-10,-5.23895e-14,-47962.8,-4.82863], Tmin=(997.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=C(O)F(172)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(1.70004,'amu*angstrom^2'), symmetry=1, barrier=(39.0874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96655,0.00185137,4.47909e-05,-7.96564e-08,4.23622e-11,-75247.4,7.80392], Tmin=(10,'K'), Tmax=(614.343,'K')), NASAPolynomial(coeffs=[2.83475,0.0173246,-1.27762e-05,4.2862e-09,-5.35304e-13,-75261.2,11.4681], Tmin=(614.343,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-625.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC(O)F""", comment="""Thermo library: CHOF_G4"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OC(O)OF(7275)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {1,S} {5,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-466.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,557,1111,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.0993805,'amu*angstrom^2'), symmetry=1, barrier=(2.28495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995613,'amu*angstrom^2'), symmetry=1, barrier=(2.28911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09967,'amu*angstrom^2'), symmetry=1, barrier=(2.29161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0312,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20014,0.0469738,-8.64586e-05,8.52523e-08,-3.18605e-11,-56012.8,17.6285], Tmin=(100,'K'), Tmax=(838.07,'K')), NASAPolynomial(coeffs=[4.1512,0.0213266,-1.13175e-05,2.22176e-09,-1.54368e-13,-55766.1,11.9828], Tmin=(838.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-466.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2sCF) + group(Cs-OsOsOsH)"""),
)

species(
    label = '[O]C(F)O[C](O)O(7309)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u1 p0 c0 {2,S} {3,S} {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-546.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,451,553,637,1069,1180,1265,1301,3056,360,370,350,276.63,1893.04,1893.47],'cm^-1')),
        HinderedRotor(inertia=(0.123921,'amu*angstrom^2'), symmetry=1, barrier=(6.69957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123381,'amu*angstrom^2'), symmetry=1, barrier=(6.70753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123872,'amu*angstrom^2'), symmetry=1, barrier=(6.70522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123172,'amu*angstrom^2'), symmetry=1, barrier=(6.70514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3635,0.0669682,-0.000119424,1.13716e-07,-4.15908e-11,-65699.4,26.4492], Tmin=(100,'K'), Tmax=(830.293,'K')), NASAPolynomial(coeffs=[5.80686,0.0272209,-1.44823e-05,2.84955e-09,-1.98496e-13,-65805.1,9.64601], Tmin=(830.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-546.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(Cs_P)"""),
)

species(
    label = 'O[C](O)O[C](O)F(7310)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-597.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3615,3650,1210,1277.5,1345,900,1000,1100,360,370,350,482,586,761,1411,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13019,0.0709689,-0.000124092,1.13271e-07,-4.01044e-11,-71787.7,27.623], Tmin=(100,'K'), Tmax=(820.763,'K')), NASAPolynomial(coeffs=[8.0007,0.02378,-1.28041e-05,2.532e-09,-1.7697e-13,-72453.9,-1.35395], Tmin=(820.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(CsFHOO) + radical(Cs_P) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = '[O]C(O)OC([O])F(7311)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {6,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-524.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,451,553,637,1069,1180,1265,1301,3056,180,1711.35,1711.98,1712.04],'cm^-1')),
        HinderedRotor(inertia=(0.202056,'amu*angstrom^2'), symmetry=1, barrier=(4.64567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201926,'amu*angstrom^2'), symmetry=1, barrier=(4.64268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201929,'amu*angstrom^2'), symmetry=1, barrier=(4.64274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87421,0.0635132,-0.000134107,1.52132e-07,-6.22113e-11,-63061.2,25.092], Tmin=(100,'K'), Tmax=(838.229,'K')), NASAPolynomial(coeffs=[-3.05336,0.0447115,-2.47376e-05,4.92216e-09,-3.44201e-13,-60748.5,56.8628], Tmin=(838.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-524.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(CsFHOO) + radical(OCOJ) + radical(O2sj(Cs-F1sO2sH))"""),
)

species(
    label = '[O]C(O)O[C](O)F(7312)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {6,S} {9,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u1 p0 c0 {1,S} {2,S} {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-575.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,482,586,761,1411,180,1232.76,1234.91,1235.92],'cm^-1')),
        HinderedRotor(inertia=(0.195817,'amu*angstrom^2'), symmetry=1, barrier=(4.50222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196378,'amu*angstrom^2'), symmetry=1, barrier=(4.51511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192517,'amu*angstrom^2'), symmetry=1, barrier=(4.42634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194958,'amu*angstrom^2'), symmetry=1, barrier=(4.48247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6363,0.0675726,-0.000139001,1.52011e-07,-6.08756e-11,-69149.3,26.2821], Tmin=(100,'K'), Tmax=(834.935,'K')), NASAPolynomial(coeffs=[-0.845718,0.0412459,-2.30446e-05,4.60103e-09,-3.22372e-13,-67402.7,45.7859], Tmin=(834.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-575.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(CsFHOO) + radical(OCOJ) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
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
    label = 'O=[C]OC(O)O(7313)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {5,S} {9,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u1 p0 c0 {1,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-555.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,1855,455,950,234.765,2056.28],'cm^-1')),
        HinderedRotor(inertia=(0.0112494,'amu*angstrom^2'), symmetry=1, barrier=(33.753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177607,'amu*angstrom^2'), symmetry=1, barrier=(6.94614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177652,'amu*angstrom^2'), symmetry=1, barrier=(6.94592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863165,'amu*angstrom^2'), symmetry=1, barrier=(33.753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0428,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99678,0.0510916,-8.62958e-05,8.4123e-08,-3.19962e-11,-66722.1,20.5218], Tmin=(100,'K'), Tmax=(812.014,'K')), NASAPolynomial(coeffs=[4.08951,0.0263316,-1.38624e-05,2.73807e-09,-1.92076e-13,-66585.6,13.795], Tmin=(812.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-555.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[O]C(=O)F(168)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-381.43,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(62.9882,'amu')),
        NonlinearRotor(inertia=([36.4335,44.7291,81.1626],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([509.533,537.43,765.261,1011.27,1202.3,1550.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[CH]O(7005)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-212.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.301364,'amu*angstrom^2'), symmetry=1, barrier=(6.92896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301175,'amu*angstrom^2'), symmetry=1, barrier=(6.92461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81665,0.0341553,-7.35611e-05,7.77e-08,-2.95104e-11,-25513.5,10.1014], Tmin=(100,'K'), Tmax=(870.402,'K')), NASAPolynomial(coeffs=[2.33151,0.0164432,-8.67075e-06,1.67642e-09,-1.14269e-13,-24673.7,16.7139], Tmin=(870.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: thermo_DFT_CCSDTF12_BAC + radical(OCJO)"""),
)

species(
    label = 'CFO(51)',
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
    label = '[O]C(O)O(7314)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {6,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-373.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.000757259,'amu*angstrom^2'), symmetry=1, barrier=(3.28267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000768368,'amu*angstrom^2'), symmetry=1, barrier=(3.34597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24526,0.0316185,-8.37877e-05,1.09317e-07,-4.73195e-11,-44885.7,14.1783], Tmin=(100,'K'), Tmax=(853.609,'K')), NASAPolynomial(coeffs=[-6.08765,0.0351523,-1.93558e-05,3.82494e-09,-2.65477e-13,-41827.8,66.3059], Tmin=(853.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-373.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + radical(OCOJ)"""),
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
    label = 'O=C(F)O[CH]O(2372)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {2,S} {3,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-595.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,482,664,788,1296,1923,265.047,558.524,2496.93],'cm^-1')),
        HinderedRotor(inertia=(0.833188,'amu*angstrom^2'), symmetry=1, barrier=(19.1566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833213,'amu*angstrom^2'), symmetry=1, barrier=(19.1572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1153,'amu*angstrom^2'), symmetry=1, barrier=(25.517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81287,0.0464932,-5.47881e-05,3.13845e-08,-6.97397e-12,-71562.3,19.9356], Tmin=(100,'K'), Tmax=(1106.43,'K')), NASAPolynomial(coeffs=[12.1838,0.00899983,-3.95791e-06,7.57392e-10,-5.37145e-14,-73857.2,-31.1472], Tmin=(1106.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(COFOO) + radical(OCJO)"""),
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
    label = '[O]C(O)OC(=O)F(7315)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {1,S} {2,S} {5,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-746.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,482,664,788,1296,1923,180,1374.69,1374.71,1374.82],'cm^-1')),
        HinderedRotor(inertia=(0.124602,'amu*angstrom^2'), symmetry=1, barrier=(2.86483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12462,'amu*angstrom^2'), symmetry=1, barrier=(2.86526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1246,'amu*angstrom^2'), symmetry=1, barrier=(2.86479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18069,0.0547526,-0.000115516,1.34376e-07,-5.64083e-11,-89786.9,23.0845], Tmin=(100,'K'), Tmax=(824.995,'K')), NASAPolynomial(coeffs=[-2.97371,0.0420119,-2.37473e-05,4.78122e-09,-3.37333e-13,-87652.4,54.742], Tmin=(824.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(COFOO) + radical(OCOJ)"""),
)

species(
    label = 'O=C(F)O[C](O)O(6084)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {2,S} {3,S} {4,S}
7 C u0 p0 c0 {1,S} {2,S} {5,D}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-769.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,482,664,788,1296,1923,180,1106.37,1107.8],'cm^-1')),
        HinderedRotor(inertia=(0.158389,'amu*angstrom^2'), symmetry=1, barrier=(3.64167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160578,'amu*angstrom^2'), symmetry=1, barrier=(3.692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160507,'amu*angstrom^2'), symmetry=1, barrier=(3.69038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158156,'amu*angstrom^2'), symmetry=1, barrier=(3.63633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68324,0.0580428,-0.000100226,9.5151e-08,-3.54502e-11,-92425.7,24.3947], Tmin=(100,'K'), Tmax=(796.084,'K')), NASAPolynomial(coeffs=[5.83287,0.0246171,-1.35493e-05,2.72248e-09,-1.928e-13,-92687.9,7.82422], Tmin=(796.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-769.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(COFOO) + radical(Cs_P)"""),
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
    E0 = (-506.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-101.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-91.0482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-196.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-107.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-158.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-89.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-201.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-171.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-174.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-142.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-164.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-447.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)OC(O)O(7264)'],
    products = ['HF(38)', 'CO2(14)', 'O=CO(1075)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.13165e+15,'s^-1'), n=-0.436116, Ea=(75.5167,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.9927400314419011, var=19.60363129417426, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'OC(O)OF(7275)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(91.1518,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)O[C](O)O(7309)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O[C](O)O[C](O)F(7310)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(O)OC([O])F(7311)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(O)O[C](O)F(7312)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'O=[C]OC(O)O(7313)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)F(168)', 'O[CH]O(7005)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CFO(51)', '[O]C(O)O(7314)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(7)', 'O=C(F)O[CH]O(2372)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]C(O)OC(=O)F(7315)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C(F)O[C](O)O(6084)'],
    products = ['O=C(F)OC(O)O(7264)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.31504e+12,'m^3/(mol*s)'), n=-2.20453, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C(F)OC(O)O(7264)'],
    products = ['O=C(O)F(172)', 'O=CO(1075)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.6407e+06,'s^-1'), n=1.67734, Ea=(134.261,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2285837836968013, var=2.96975919141528, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #2820',
    isomers = [
        'O=C(F)OC(O)O(7264)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'O=CO(1075)'),
        ('O=C(O)F(172)', 'O=CO(1075)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2820',
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

