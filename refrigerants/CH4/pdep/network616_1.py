species(
    label = '[CH]=CO[C](O)O(1146)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-77.2972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,317.18,318.728],'cm^-1')),
        HinderedRotor(inertia=(0.0016867,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195889,'amu*angstrom^2'), symmetry=1, barrier=(13.8157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19662,'amu*angstrom^2'), symmetry=1, barrier=(13.8131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193456,'amu*angstrom^2'), symmetry=1, barrier=(13.8188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16452,0.0595778,-7.19529e-05,4.18451e-08,-9.36482e-12,-9192.16,23.1347], Tmin=(100,'K'), Tmax=(1103.79,'K')), NASAPolynomial(coeffs=[14.9903,0.00947503,-3.86589e-06,7.22191e-10,-5.08388e-14,-12244.3,-44.9322], Tmin=(1103.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.2972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = 'O=C(O)O(747)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,421.686,423.663,424.543,424.784],'cm^-1')),
        HinderedRotor(inertia=(0.00093875,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344145,'amu*angstrom^2'), symmetry=1, barrier=(43.4791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.87,'J/mol'), sigma=(5.42485,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.51 K, Pc=51.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98124,0.0114379,2.77981e-05,-4.94746e-08,2.0802e-11,-75136.3,10.9454], Tmin=(100,'K'), Tmax=(956.736,'K')), NASAPolynomial(coeffs=[11.8951,0.00171649,-1.48201e-07,9.26559e-11,-1.38667e-14,-78102.6,-38.253], Tmin=(956.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(O2s-(Cds-O2d)H) + group(Cds-OdOsOs)"""),
)

species(
    label = 'C2H2(22)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O[C]O(630)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {5,S}
3 C u2 p0 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (50.0085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.863664,'amu*angstrom^2'), symmetry=1, barrier=(19.8573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.864861,'amu*angstrom^2'), symmetry=1, barrier=(19.8849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68163,0.0133808,3.05183e-05,-7.11677e-08,3.49255e-11,6076.95,10.6719], Tmin=(100,'K'), Tmax=(886.889,'K')), NASAPolynomial(coeffs=[19.125,-0.0179815,1.11743e-05,-2.21414e-09,1.50426e-13,1477.02,-76.1735], Tmin=(886.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.0085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(99.7737,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C[O](171)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u1 p0 c0 {2,D} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (150.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32152,0.00995633,8.2355e-06,-1.92995e-08,8.375e-12,18150.6,9.47658], Tmin=(100,'K'), Tmax=(969.359,'K')), NASAPolynomial(coeffs=[7.73862,0.00393429,-1.33187e-06,2.69027e-10,-2.15556e-14,16720.8,-14.6537], Tmin=(969.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'OC1(O)C=CO1(1148)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {4,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {1,S} {5,D} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-382.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71997,0.0338958,1.76929e-05,-6.31339e-08,3.14421e-11,-45904.6,17.5015], Tmin=(100,'K'), Tmax=(910.379,'K')), NASAPolynomial(coeffs=[19.3437,-0.000493393,3.43013e-06,-7.51523e-10,4.85754e-14,-50897.2,-75.665], Tmin=(910.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-382.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=COC(=O)O(1147)',
    structure = adjacencyList("""1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {6,D} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-531.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4433,0.041636,-3.00299e-06,-3.74019e-08,2.00431e-11,-63829.3,19.2274], Tmin=(100,'K'), Tmax=(956.19,'K')), NASAPolynomial(coeffs=[18.8941,0.00442633,-7.78836e-07,1.94166e-10,-2.16617e-14,-68802.8,-72.7369], Tmin=(956.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-531.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=[CH](149)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,D} {3,S}
2 C u1 p0 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1227.31,2941.71],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88058,-0.000843434,2.04094e-05,-2.51944e-08,9.47575e-12,71059.3,4.72115], Tmin=(100,'K'), Tmax=(935.367,'K')), NASAPolynomial(coeffs=[4.92141,0.00358272,-9.24433e-07,1.57271e-10,-1.19539e-14,70476.2,-2.30656], Tmin=(935.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=COC(=O)O(1149)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 C u1 p0 c0 {4,D} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-284.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,374.999,375.016,375.04,375.042,375.073,375.1],'cm^-1')),
        HinderedRotor(inertia=(0.263117,'amu*angstrom^2'), symmetry=1, barrier=(26.2849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263036,'amu*angstrom^2'), symmetry=1, barrier=(26.2791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263275,'amu*angstrom^2'), symmetry=1, barrier=(26.2794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37047,0.0451482,-1.96271e-05,-2.00565e-08,1.41396e-11,-34109.6,19.8919], Tmin=(100,'K'), Tmax=(954.202,'K')), NASAPolynomial(coeffs=[19.1316,0.00159797,2.52701e-07,-4.03358e-12,-6.84017e-15,-38906,-72.3346], Tmin=(954.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[O][C](O)O(746)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-168.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,2095.66],'cm^-1')),
        HinderedRotor(inertia=(0.00195616,'amu*angstrom^2'), symmetry=1, barrier=(6.09647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00195478,'amu*angstrom^2'), symmetry=1, barrier=(6.09289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27975,0.0323757,-9.4266e-05,1.21485e-07,-5.16465e-11,-20202.9,14.7429], Tmin=(100,'K'), Tmax=(861.003,'K')), NASAPolynomial(coeffs=[-6.25243,0.0327858,-1.85451e-05,3.67222e-09,-2.54021e-13,-16935.2,68.7475], Tmin=(861.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + radical(OCOJ) + radical(Cs_P)"""),
)

species(
    label = 'C#CO[C](O)O(1150)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {4,S} {8,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 C u0 p0 c0 {1,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-115.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,360,370,350,2175,525,750,770,3400,2100,305.605,306.47],'cm^-1')),
        HinderedRotor(inertia=(0.00180387,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26601,'amu*angstrom^2'), symmetry=1, barrier=(17.4508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262596,'amu*angstrom^2'), symmetry=1, barrier=(17.4234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265774,'amu*angstrom^2'), symmetry=1, barrier=(17.4394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40883,0.059579,-8.35502e-05,5.74192e-08,-1.53956e-11,-13748.2,19.319], Tmin=(100,'K'), Tmax=(917.444,'K')), NASAPolynomial(coeffs=[12.2544,0.0122947,-6.24445e-06,1.24649e-09,-8.93351e-14,-15738.3,-32.0705], Tmin=(917.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cs_P)"""),
)

species(
    label = '[CH]=COC([O])O(1151)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
5  C u0 p0 c0 {1,S} {6,D} {8,S}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-55.1457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,1315.93],'cm^-1')),
        HinderedRotor(inertia=(0.140243,'amu*angstrom^2'), symmetry=1, barrier=(3.22445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136697,'amu*angstrom^2'), symmetry=1, barrier=(3.14293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17175,'amu*angstrom^2'), symmetry=1, barrier=(26.9409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74947,0.055253,-8.37772e-05,7.72193e-08,-2.92099e-11,-6557.02,21.5122], Tmin=(100,'K'), Tmax=(757.12,'K')), NASAPolynomial(coeffs=[5.38538,0.0282608,-1.48807e-05,2.97649e-09,-2.11759e-13,-6884.51,6.45583], Tmin=(757.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.1457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]O[C](O)O(1152)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {4,S} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  C u0 p0 c0 {6,D} {7,S} {8,S}
6  C u1 p0 c0 {1,S} {5,D}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-84.6492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61677,0.0557041,-7.42803e-05,5.23579e-08,-1.47817e-11,-10098,24.2457], Tmin=(100,'K'), Tmax=(864.131,'K')), NASAPolynomial(coeffs=[9.79362,0.0178547,-8.58057e-06,1.67217e-09,-1.1814e-13,-11511.2,-14.009], Tmin=(864.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.6492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OC(O)O(1153)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {4,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
5  C u1 p0 c0 {1,S} {6,D}
6  C u1 p0 c0 {5,D} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-42.7994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,1685,370,3120,650,792.5,1650,340.62,340.67],'cm^-1')),
        HinderedRotor(inertia=(0.00145329,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11703,'amu*angstrom^2'), symmetry=1, barrier=(9.64123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117044,'amu*angstrom^2'), symmetry=1, barrier=(9.64114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117057,'amu*angstrom^2'), symmetry=1, barrier=(9.64012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47718,0.058846,-8.17929e-05,5.92886e-08,-1.70778e-11,-5059.73,24.4608], Tmin=(100,'K'), Tmax=(850.141,'K')), NASAPolynomial(coeffs=[10.3412,0.0171421,-8.21368e-06,1.59194e-09,-1.11876e-13,-6566.93,-16.8642], Tmin=(850.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.7994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]OC(=O)O(1154)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u1 p0 c0 {1,S} {6,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u1 p0 c0 {4,S} {8,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
"""),
    E0 = (-239.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13143,0.0552635,-5.11696e-05,1.88027e-08,-1.4021e-12,-28693.2,23.548], Tmin=(100,'K'), Tmax=(1035.05,'K')), NASAPolynomial(coeffs=[16.5353,0.00912162,-3.70111e-06,7.24447e-10,-5.37077e-14,-32599.1,-54.7615], Tmin=(1035.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsOs) + radical(CCsJOC(O)) + radical(CJCO)"""),
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
    E0 = (-52.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (225.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-44.5237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-18.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (26.5616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (37.0225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (87.7562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (128.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (447.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (127.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (52.6288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (135.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-19.7678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CO[C](O)O(1146)'],
    products = ['O=C(O)O(747)', 'C2H2(22)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O[C]O(630)', '[CH]=C[O](171)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CO[C](O)O(1146)'],
    products = ['OC1(O)C=CO1(1148)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_NDMustO;CdsinglepriH_rad_out]
Euclidian distance = 3.7416573867739413
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CO[C](O)O(1146)'],
    products = ['C=COC(=O)O(1147)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(O)O(747)', '[CH]=[CH](149)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.049788,'m^3/(mol*s)'), n=2.39813, Ea=(36.4204,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3034033367420285, var=0.7726965170871078, Tref=1000.0, N=21, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(6)', '[CH]=COC(=O)O(1149)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(18.3,'m^3/(mol*s)'), n=1.99, Ea=(85.215,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-2C-R_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_1COS->O_2CS->C_Ext-2C-R_Ext-2C-R_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C](O)O(746)', 'C2H2(22)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(13.5425,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(6)', 'C#CO[C](O)O(1150)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1108.07,'m^3/(mol*s)'), n=1.405, Ea=(7.71778,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C](O)O(746)', '[CH]=[CH](149)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=COC([O])O(1151)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CO[C](O)O(1146)'],
    products = ['C=[C]O[C](O)O(1152)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OC(O)O(1153)'],
    products = ['[CH]=CO[C](O)O(1146)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3e+08,'s^-1'), n=1.23, Ea=(154.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_NDMustO] for rate rule [R3H_SS_O;Cd_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CO[C](O)O(1146)'],
    products = ['[CH2][CH]OC(=O)O(1154)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #616',
    isomers = [
        '[CH]=CO[C](O)O(1146)',
    ],
    reactants = [
        ('O=C(O)O(747)', 'C2H2(22)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #616',
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

