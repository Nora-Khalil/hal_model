species(
    label = '[CH2]C(F)C1(F)[CH]C1(F)F(11631)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {11,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-351.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,2950,1000,3000,3100,440,815,1455,1000,180,180,1590.22],'cm^-1')),
        HinderedRotor(inertia=(0.0270856,'amu*angstrom^2'), symmetry=1, barrier=(48.4586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156111,'amu*angstrom^2'), symmetry=1, barrier=(3.5893,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.590429,0.0812261,-0.000101263,6.93294e-08,-1.95424e-11,-42197.8,27.2484], Tmin=(100,'K'), Tmax=(855.071,'K')), NASAPolynomial(coeffs=[10.9122,0.0329423,-1.65631e-05,3.29342e-09,-2.35614e-13,-43963,-20.9321], Tmin=(855.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-351.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'CH2CHF(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-153.05,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0219,'amu')),
        NonlinearRotor(inertia=([7.59478,47.6085,55.2033],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([482.184,740.114,880.476,949.569,983.363,1189.2,1343.65,1421.6,1725.69,3171.53,3191.61,3269.97],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2263.2,'J/mol'), sigma=(4.322,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09164,-0.0073724,7.45741e-05,-1.12982e-07,5.61696e-11,-18407.3,6.78145], Tmin=(10,'K'), Tmax=(619.705,'K')), NASAPolynomial(coeffs=[1.44203,0.0189088,-1.12569e-05,3.25441e-09,-3.64262e-13,-18255.1,16.8744], Tmin=(619.705,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-153.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CC1(F)F(974)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,D} {7,S}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-329.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,323,467,575,827,1418,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2932.45,'J/mol'), sigma=(5.0155,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.04 K, Pc=52.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90533,0.00560985,9.08669e-05,-1.89954e-07,1.15706e-10,-39564.3,10.133], Tmin=(10,'K'), Tmax=(565.045,'K')), NASAPolynomial(coeffs=[4.14369,0.0253568,-1.84551e-05,6.16339e-09,-7.67829e-13,-39933.5,6.09131], Tmin=(565.045,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-329.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C(F)C1[C](F)C1(F)F(11632)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {9,S} {11,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-352.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,212,367,445,1450,3000,3100,440,815,1455,1000,180,646.533,646.583,646.601,646.613],'cm^-1')),
        HinderedRotor(inertia=(0.00854557,'amu*angstrom^2'), symmetry=1, barrier=(2.53458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00853864,'amu*angstrom^2'), symmetry=1, barrier=(2.53324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.5,'J/mol'), sigma=(5.6033,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.45 K, Pc=40.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813301,0.0731824,-7.61387e-05,4.1789e-08,-9.40441e-12,-42337.8,26.69], Tmin=(100,'K'), Tmax=(1059.07,'K')), NASAPolynomial(coeffs=[12.3944,0.0294427,-1.41899e-05,2.79413e-09,-1.99635e-13,-44790.9,-29.8471], Tmin=(1059.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-352.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsCsFH) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]CC1(F)[CH]C1(F)F(11633)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-350.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,2950,1000,334,575,1197,1424,3202,180,966.502,968.007,972.44],'cm^-1')),
        HinderedRotor(inertia=(0.200268,'amu*angstrom^2'), symmetry=1, barrier=(4.60456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200744,'amu*angstrom^2'), symmetry=1, barrier=(4.6155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3179.54,'J/mol'), sigma=(5.92973,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=496.64 K, Pc=34.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.551187,0.082478,-0.000109286,8.12812e-08,-2.49322e-11,-42018.1,27.9766], Tmin=(100,'K'), Tmax=(788.593,'K')), NASAPolynomial(coeffs=[10.0989,0.0340481,-1.71649e-05,3.40133e-09,-2.42217e-13,-43523.9,-15.8176], Tmin=(788.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-350.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Csj(Cs-CsHH)(F1s)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'F[CH]C1(F)[CH]C1(F)F(11484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-317.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,2950,1000,334,575,1197,1424,3202,541.842,542.001,542.237],'cm^-1')),
        HinderedRotor(inertia=(0.000573644,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15803,0.0639692,-7.44679e-05,4.36093e-08,-1.01166e-11,-38049.1,23.2582], Tmin=(100,'K'), Tmax=(1048.07,'K')), NASAPolynomial(coeffs=[13.2155,0.0179519,-8.60838e-06,1.71718e-09,-1.24049e-13,-40576.5,-35.4784], Tmin=(1048.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Csj(Cs-F1sCsCs)(F1s)(H)_1901_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1CC2C(F)(F)C12F(11735)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {6,S} {8,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-679.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09353,0.0618585,-4.84493e-05,1.83108e-08,-2.79017e-12,-81605.3,19.168], Tmin=(100,'K'), Tmax=(1519.68,'K')), NASAPolynomial(coeffs=[15.089,0.0250209,-1.20892e-05,2.36027e-09,-1.66201e-13,-85859.1,-54.2095], Tmin=(1519.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'C=C(F)C1(F)CC1(F)F(11800)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {8,D} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-679.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961285,0.069389,-6.82334e-05,3.49789e-08,-7.36456e-12,-81567.5,23.6622], Tmin=(100,'K'), Tmax=(1124.3,'K')), NASAPolynomial(coeffs=[12.4165,0.0286352,-1.38629e-05,2.74032e-09,-1.96195e-13,-84143.4,-32.9451], Tmin=(1124.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CdCsCdF) + group(Cds-CdsHH) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = '[CH2]C(F)[C](F)C=C(F)F(10532)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u1 p0 c0 {2,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u1 p0 c0 {5,S} {12,S} {13,S}
9  C u0 p0 c0 {3,S} {4,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-486.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,346,659,817,1284,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413,180,592.197],'cm^-1')),
        HinderedRotor(inertia=(0.773871,'amu*angstrom^2'), symmetry=1, barrier=(17.7928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.32915,'amu*angstrom^2'), symmetry=1, barrier=(53.5517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219491,'amu*angstrom^2'), symmetry=1, barrier=(5.04653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582544,0.0813724,-0.000100951,6.80347e-08,-1.88284e-11,-58357.2,28.0076], Tmin=(100,'K'), Tmax=(870.87,'K')), NASAPolynomial(coeffs=[11.2781,0.0322467,-1.63368e-05,3.26085e-09,-2.3387e-13,-60220.1,-22.1135], Tmin=(870.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-486.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCFF) + radical(CsCdCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
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
    label = '[CH2]C(F)C1=CC1(F)F(11801)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {4,S}
4  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {10,S}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-216.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,174,267,591,721,1107,1278,1348,3273,2950,1000,3000,3100,440,815,1455,1000,360.147,360.156,360.162,360.182,360.184],'cm^-1')),
        HinderedRotor(inertia=(0.31059,'amu*angstrom^2'), symmetry=1, barrier=(28.5925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310605,'amu*angstrom^2'), symmetry=1, barrier=(28.5922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23528,0.0624499,-5.94206e-05,2.87068e-08,-5.66033e-12,-25996.4,23.552], Tmin=(100,'K'), Tmax=(1197.06,'K')), NASAPolynomial(coeffs=[12.4974,0.0248174,-1.22647e-05,2.44482e-09,-1.75683e-13,-28692.7,-32.8071], Tmin=(1197.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCCFF) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd-Cd-Cs(F)(F)) + radical(Csj(Cs-F1sCdH)(H)(H))"""),
)

species(
    label = '[CH2][CH]F(1279)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (116.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,334,575,1197,1424,3202],'cm^-1')),
        HinderedRotor(inertia=(0.0021411,'amu*angstrom^2'), symmetry=1, barrier=(6.24533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56892,0.0121147,-4.46307e-06,3.52998e-10,6.36808e-14,13988.6,11.2355], Tmin=(100,'K'), Tmax=(2106.76,'K')), NASAPolynomial(coeffs=[7.95456,0.00674556,-2.74612e-06,4.76062e-10,-2.99992e-14,11484.4,-14.7484], Tmin=(2106.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sH)"""),
)

species(
    label = '[CH2]C(F)C1(F)C=C1F(11802)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {10,S}
8  C u1 p0 c0 {5,S} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-165.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,259,529,569,1128,1321,1390,3140,323,467,575,827,1418,2950,1000,3000,3100,440,815,1455,1000,180,1404.76],'cm^-1')),
        HinderedRotor(inertia=(0.449215,'amu*angstrom^2'), symmetry=1, barrier=(10.3283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17419,'amu*angstrom^2'), symmetry=1, barrier=(26.997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902805,0.074343,-0.0001048,8.50488e-08,-2.85711e-11,-19814.8,23.416], Tmin=(100,'K'), Tmax=(721.692,'K')), NASAPolynomial(coeffs=[8.59736,0.031704,-1.61942e-05,3.21455e-09,-2.28512e-13,-20925.6,-11.1976], Tmin=(721.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cd(F)-Cd) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'C=CC1(F)[CH]C1(F)F(11803)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
6  C u1 p0 c0 {4,S} {5,S} {9,S}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-244.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,317.043,317.273,317.42],'cm^-1')),
        HinderedRotor(inertia=(0.00167715,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73339,0.0544628,-4.71098e-05,2.20859e-08,-4.47193e-12,-29349.5,23.1415], Tmin=(100,'K'), Tmax=(1122.79,'K')), NASAPolynomial(coeffs=[8.50467,0.0303396,-1.48821e-05,2.95031e-09,-2.11202e-13,-30870,-10.3102], Tmin=(1122.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C1(F)F(11362)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-103.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2950,1000,212,367,445,1450,2726.48,2726.56],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00781,0.0264579,-1.29283e-05,-2.29098e-08,2.96552e-11,-12413.9,18.9514], Tmin=(100,'K'), Tmax=(451.813,'K')), NASAPolynomial(coeffs=[4.42599,0.0197907,-1.03421e-05,2.11919e-09,-1.54582e-13,-12602.1,12.5711], Tmin=(451.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
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
    label = 'C=C(F)C1(F)[CH]C1(F)F(11395)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u1 p0 c0 {5,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {8,D} {11,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-453.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,2950,1000,323,467,575,827,1418,2950,3100,1380,975,1025,1650,180,180,1244.15],'cm^-1')),
        HinderedRotor(inertia=(0.236212,'amu*angstrom^2'), symmetry=1, barrier=(5.43097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09285,0.0701706,-8.63729e-05,6.03447e-08,-1.76828e-11,-54404.6,25.5513], Tmin=(100,'K'), Tmax=(817.586,'K')), NASAPolynomial(coeffs=[8.88774,0.0320339,-1.64037e-05,3.29038e-09,-2.36552e-13,-55679.2,-10.4845], Tmin=(817.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-453.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CdCsCdF) + group(Cds-CdsHH) + ring(Cs-Cs(F)(F)-Cs) + radical(cyclopropane) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = '[CH2]C(F)[C]1C=C1F(11804)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {4,S} {5,D} {9,S}
7  C u1 p0 c0 {3,S} {10,S} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (175.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,271,519,563,612,1379,2950,1000,3000,3100,440,815,1455,1000,180,1021.86,1022.01,1022.01,1022.15],'cm^-1')),
        HinderedRotor(inertia=(0.199484,'amu*angstrom^2'), symmetry=1, barrier=(4.58653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199399,'amu*angstrom^2'), symmetry=1, barrier=(4.58457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63898,0.0560734,-6.51154e-05,4.53575e-08,-1.34439e-11,21134.1,18.1799], Tmin=(100,'K'), Tmax=(806.124,'K')), NASAPolynomial(coeffs=[7.23867,0.0282874,-1.34123e-05,2.59862e-09,-1.83181e-13,20231.3,-7.62851], Tmin=(806.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs-Cd(F)-Cd) + radical(Allyl_T) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = '[CH2]C=C1[CH]C1(F)F(11805)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  F u0 p3 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u1 p0 c0 {3,S} {4,S} {8,S}
6  C u0 p0 c0 {4,D} {7,S} {9,S}
7  C u1 p0 c0 {6,S} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (56.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2950,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,647.41,647.466,647.826,647.942,647.971],'cm^-1')),
        HinderedRotor(inertia=(0.000402209,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7323,0.0374806,5.89306e-06,-3.76851e-08,1.77254e-11,6831.35,20.3391], Tmin=(100,'K'), Tmax=(978.103,'K')), NASAPolynomial(coeffs=[14.3494,0.0154438,-5.64656e-06,1.08003e-09,-8.06946e-14,2949.15,-47.4801], Tmin=(978.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(Allyl_S) + radical(Allyl_P)"""),
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
    label = '[CH2]C(F)=C1[CH]C1(F)F(11397)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  C u1 p0 c0 {4,S} {5,S} {9,S}
7  C u0 p0 c0 {3,S} {5,D} {8,S}
8  C u1 p0 c0 {7,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-158.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,2950,1000,271,519,563,612,1379,3000,3100,440,815,1455,1000,675.879,675.887,675.896,675.918,675.931],'cm^-1')),
        HinderedRotor(inertia=(0.000368961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06977,0.055723,-4.1325e-05,9.66043e-09,1.09923e-12,-18894.5,23.106], Tmin=(100,'K'), Tmax=(1058.23,'K')), NASAPolynomial(coeffs=[15.484,0.0162711,-6.71106e-06,1.27769e-09,-9.15979e-14,-22786.9,-51.2268], Tmin=(1058.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cs-Cs(F)(F)-Cd(Cd)) + radical(Allyl_S) + radical(Csj(Cd-F1sCd)(H)(H))"""),
)

species(
    label = 'C[C](F)C1(F)[CH]C1(F)F(11806)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {9,S} {10,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {6,S} {13,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-368.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,1000,212,367,445,1450,180,1375.68,1378.42],'cm^-1')),
        HinderedRotor(inertia=(0.177771,'amu*angstrom^2'), symmetry=1, barrier=(4.0873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180022,'amu*angstrom^2'), symmetry=1, barrier=(4.13906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935473,0.0773173,-8.82547e-05,3.6666e-08,8.51372e-12,-44169.8,26.6821], Tmin=(100,'K'), Tmax=(546.02,'K')), NASAPolynomial(coeffs=[8.64844,0.036368,-1.84899e-05,3.65657e-09,-2.59241e-13,-45243.9,-7.98463], Tmin=(546.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2][C](F)C1(F)CC1(F)F(11807)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u1 p0 c0 {4,S} {5,S} {9,S}
9  C u1 p0 c0 {8,S} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-395.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.97626,0.0760133,-8.6329e-05,3.76374e-08,6.15561e-12,-47511.4,27.5281], Tmin=(100,'K'), Tmax=(551.315,'K')), NASAPolynomial(coeffs=[8.45213,0.0364077,-1.8389e-05,3.63102e-09,-2.57393e-13,-48558.1,-6.10428], Tmin=(551.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-395.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2]C(F)[C]1C(F)C1(F)F(10557)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-392.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59622,0.0653702,-3.26346e-05,-9.19475e-08,1.16929e-10,-47093.2,25.3294], Tmin=(100,'K'), Tmax=(451.77,'K')), NASAPolynomial(coeffs=[7.42747,0.0374521,-1.86703e-05,3.62769e-09,-2.53042e-13,-47862.1,-0.84772], Tmin=(451.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-392.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FCC(F)[C]1[CH]C1(F)F(11808)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {6,S} {9,S}
9  C u1 p0 c0 {6,S} {8,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-359.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997247,0.0745485,-8.62692e-05,4.66949e-08,-3.45957e-12,-43145.4,27.361], Tmin=(100,'K'), Tmax=(573.688,'K')), NASAPolynomial(coeffs=[8.08325,0.0368893,-1.85191e-05,3.65832e-09,-2.59866e-13,-44151.8,-4.57227], Tmin=(573.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-359.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(H)_ring)"""),
)

species(
    label = '[CH2]C(F)C1(F)[C](F)C1F(11809)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-327.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,207,311,436,622,447,691,1074,1182,1276,1366,1311,1469,3045,3235,212,367,445,1450,3000,3100,440,815,1455,1000,365.547],'cm^-1')),
        HinderedRotor(inertia=(0.00125969,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00125964,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.717188,0.0767254,-8.84895e-05,5.53134e-08,-1.4183e-11,-39235.7,28.2568], Tmin=(100,'K'), Tmax=(937.566,'K')), NASAPolynomial(coeffs=[11.5165,0.030651,-1.47748e-05,2.89721e-09,-2.06169e-13,-41260.7,-23.1475], Tmin=(937.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sH)(F1s)_ring) + radical(Csj(Cs-F1sCsH)(H)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FCC(F)C1(F)[CH][C]1F(11810)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {10,S}
7  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
8  C u1 p0 c0 {4,S} {5,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-309.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,250,417,511,1155,1315,1456,3119,528,1116,1182,1331,1402,1494,3075,3110,212,367,445,1450,2950,1000,436.568,436.593,436.882,436.934],'cm^-1')),
        HinderedRotor(inertia=(0.000885603,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350051,'amu*angstrom^2'), symmetry=1, barrier=(47.2972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728849,0.0776142,-8.83987e-05,5.43252e-08,-1.38004e-11,-37127.4,26.7507], Tmin=(100,'K'), Tmax=(940.986,'K')), NASAPolynomial(coeffs=[11.3639,0.0324055,-1.63316e-05,3.26657e-09,-2.35028e-13,-39128.8,-23.9102], Tmin=(940.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-F1sCsCs)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sH)(H)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH2][CH]C1(F)C(F)C1(F)F(11811)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
7  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
8  C u1 p0 c0 {5,S} {9,S} {11,S}
9  C u1 p0 c0 {8,S} {12,S} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-338.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,222,329,445,522,589,1214,1475,250,417,511,1155,1315,1456,3119,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1392.26],'cm^-1')),
        HinderedRotor(inertia=(0.00183616,'amu*angstrom^2'), symmetry=1, barrier=(2.5257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00183615,'amu*angstrom^2'), symmetry=1, barrier=(2.52567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14454,0.068641,-7.59331e-05,4.91878e-08,-1.36593e-11,-40594.2,29.2061], Tmin=(100,'K'), Tmax=(852.432,'K')), NASAPolynomial(coeffs=[8.22679,0.0354083,-1.74553e-05,3.4544e-09,-2.46893e-13,-41801.7,-3.83101], Tmin=(852.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-338.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Cs_S) + radical(RCCJ) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC[CH]C1(F)[CH]C1(F)F(11750)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {13,S}
9  C u1 p0 c0 {5,S} {7,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-347.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863006,0.0759211,-9.45935e-05,6.87697e-08,-2.11749e-11,-41736.6,28.2772], Tmin=(100,'K'), Tmax=(777.428,'K')), NASAPolynomial(coeffs=[8.49737,0.0366387,-1.87958e-05,3.76706e-09,-2.70489e-13,-42923.5,-6.63167], Tmin=(777.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Cs_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-165.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-8.66447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-7.3561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (250.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-157.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-102.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-74.4067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (98.5944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-18.1627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (136.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (79.4907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-61.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-47.0107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (198.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (359.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (253.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (55.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (2.40003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-8.24888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-3.62888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (57.4049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (35.1839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (97.8191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (69.3873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (23.8265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['CH2CHF(56)', 'FC1=CC1(F)F(974)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['[CH2]C(F)C1[C](F)C1(F)F(11632)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['F[CH]CC1(F)[CH]C1(F)F(11633)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(18)', 'F[CH]C1(F)[CH]C1(F)F(11484)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['FC1CC2C(F)(F)C12F(11735)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['C=C(F)C1(F)CC1(F)F(11800)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)[C](F)C=C(F)F(10532)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', '[CH2]C(F)C1=CC1(F)F(11801)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(56.8192,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]F(1279)', 'FC1=CC1(F)F(974)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.47543e-05,'m^3/(mol*s)'), n=2.68942, Ea=(8.79827,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5315216239348307, var=1.6684324977447595, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-1R!H-R_Ext-3C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing_Ext-1R!H-R_Ext-3C-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2]C(F)C1(F)C=C1F(11802)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(43.1084,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C=CC1(F)[CH]C1(F)F(11803)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(65.4252,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2CHF(56)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(8.69171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'C=C(F)C1(F)[CH]C1(F)F(11395)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0317097,'m^3/(mol*s)'), n=2.69193, Ea=(8.50917,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.02638351801464579, var=1.2041797846323272, Tref=1000.0, N=199, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]F(1279)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F2(78)', '[CH2]C(F)[C]1C=C1F(11804)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(7.16775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(78)', '[CH2]C=C1[CH]C1(F)F(11805)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(20.0197,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', '[CH2]C(F)=C1[CH]C1(F)F(11397)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(309.139,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[C](F)C1(F)[CH]C1(F)F(11806)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.465e+11,'s^-1'), n=0, Ea=(184.64,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 440 used for R2H_S;C_rad_out_noH;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['[CH2][C](F)C1(F)CC1(F)F(11807)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.2544e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['[CH2]C(F)[C]1C(F)C1(F)F(10557)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(162.354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['FCC(F)[C]1[CH]C1(F)F(11808)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(223.388,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(F)C1(F)[C](F)C1F(11809)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(176.514,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['FCC(F)C1(F)[CH][C]1F(11810)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000550858,'s^-1'), n=4.50663, Ea=(221.608,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C1(F)C(F)C1(F)F(11811)'],
    products = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(221.874,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(F)C1(F)[CH]C1(F)F(11631)'],
    products = ['FC[CH]C1(F)[CH]C1(F)F(11750)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(189.809,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #3460',
    isomers = [
        '[CH2]C(F)C1(F)[CH]C1(F)F(11631)',
    ],
    reactants = [
        ('CH2CHF(56)', 'FC1=CC1(F)F(974)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3460',
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

