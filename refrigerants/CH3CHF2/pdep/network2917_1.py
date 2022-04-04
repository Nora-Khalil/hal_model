species(
    label = '[CH]=CC(F)[C]=C(7411)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {6,D} {8,S}
4  C u0 p0 c0 {5,D} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {3,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (371.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.972207,'amu*angstrom^2'), symmetry=1, barrier=(22.353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970968,'amu*angstrom^2'), symmetry=1, barrier=(22.3245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15494,0.0475737,-3.7017e-05,-1.00547e-08,2.7552e-11,44698.8,19.4013], Tmin=(100,'K'), Tmax=(502.947,'K')), NASAPolynomial(coeffs=[5.84489,0.0290944,-1.43147e-05,2.81433e-09,-1.99613e-13,44190.2,2.76898], Tmin=(502.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(23)',
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
    label = 'C=C=CF(6101)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (9.45959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1537.31],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2831,'J/mol'), sigma=(4.74838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=442.20 K, Pc=60 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96414,0.0021816,6.68618e-05,-1.27032e-07,7.57672e-11,1138.94,8.06307], Tmin=(10,'K'), Tmax=(518.59,'K')), NASAPolynomial(coeffs=[2.17835,0.0231528,-1.46137e-05,4.46872e-09,-5.27221e-13,1227.38,14.5728], Tmin=(518.59,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(9.45959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CDCDCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=CC([CH2])=CF(7410)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,D}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {2,D} {8,S}
6  C u1 p0 c0 {3,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (263.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35286,0.0492452,-2.83735e-05,-5.69945e-09,7.79623e-12,31751.3,19.6352], Tmin=(100,'K'), Tmax=(953.529,'K')), NASAPolynomial(coeffs=[15.3022,0.0121566,-3.73811e-06,6.44391e-10,-4.61289e-14,28116.9,-52.1069], Tmin=(953.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[C]=C(395)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u2 p0 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94832,-0.00181671,2.07132e-05,-2.367e-08,8.455e-12,72206.8,5.31544], Tmin=(100,'K'), Tmax=(953.655,'K')), NASAPolynomial(coeffs=[4.20269,0.00473429,-1.57308e-06,2.85913e-10,-2.08649e-14,71811.9,2.28411], Tmin=(953.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]C=CF(270)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 C u2 p0 c0 {2,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (184.694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,422.185,423.738,440.479],'cm^-1')),
        HinderedRotor(inertia=(0.393389,'amu*angstrom^2'), symmetry=1, barrier=(51.4296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2822.22,'J/mol'), sigma=(4.85549,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=440.82 K, Pc=55.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03729,0.0164789,1.1221e-05,-1.97621e-08,6.93534e-12,22252.1,13.2192], Tmin=(100,'K'), Tmax=(1095.7,'K')), NASAPolynomial(coeffs=[5.52395,0.018931,-7.92013e-06,1.48778e-09,-1.04334e-13,21015.1,-2.16317], Tmin=(1095.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1C=CC1F(7714)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {6,D}
4  C u0 p0 c0 {2,S} {5,D} {8,S}
5  C u0 p0 c0 {3,S} {4,D} {9,S}
6  C u0 p0 c0 {3,D} {10,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (17.1894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05913,0.029654,2.06686e-05,-4.99249e-08,2.14389e-11,2148.96,13.082], Tmin=(100,'K'), Tmax=(977.683,'K')), NASAPolynomial(coeffs=[13.3932,0.0149868,-5.46929e-06,1.06558e-09,-8.08957e-14,-1582.48,-49.0914], Tmin=(977.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.1894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C=C(F)C=C(7715)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,D} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {6,D}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (51.0767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71758,0.0470898,-3.90732e-05,1.68702e-08,-2.9606e-12,6227.77,18.0137], Tmin=(100,'K'), Tmax=(1346.83,'K')), NASAPolynomial(coeffs=[11.1291,0.0191383,-7.94281e-06,1.46098e-09,-1.00325e-13,3692.64,-30.1935], Tmin=(1346.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.0767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCddCF) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH]=CC=C=C(6456)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {4,D} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {6,S}
3  C u0 p0 c0 {4,D} {8,S} {9,S}
4  C u0 p0 c0 {1,D} {3,D}
5  C u1 p0 c0 {2,D} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (481.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.29416,'amu*angstrom^2'), symmetry=1, barrier=(29.7552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0931,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98895,0.0358189,-9.42096e-06,-1.77109e-08,1.06169e-11,58002.7,16.4759], Tmin=(100,'K'), Tmax=(954.569,'K')), NASAPolynomial(coeffs=[12.5831,0.0117991,-3.69124e-06,6.46389e-10,-4.65899e-14,55051.9,-39.0039], Tmin=(954.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[CH](264)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88058,-0.000843366,2.04091e-05,-2.5194e-08,9.47559e-12,71059.3,4.72117], Tmin=(100,'K'), Tmax=(935.371,'K')), NASAPolynomial(coeffs=[4.92143,0.00358269,-9.24415e-07,1.57267e-10,-1.19536e-14,70476.2,-2.30666], Tmin=(935.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=CC(F)=C=C(7716)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {5,D}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u0 p0 c0 {5,D} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {4,D}
6  C u1 p0 c0 {3,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (298.173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.41832,'amu*angstrom^2'), symmetry=1, barrier=(32.6099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77963,0.0490009,-5.01199e-05,2.71344e-08,-5.94796e-12,35941.8,18.196], Tmin=(100,'K'), Tmax=(1096.55,'K')), NASAPolynomial(coeffs=[10.2263,0.018189,-7.97121e-06,1.5092e-09,-1.057e-13,34089.3,-23.3329], Tmin=(1096.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCddCF) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH]F(6102)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (196.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.900552,'amu*angstrom^2'), symmetry=1, barrier=(20.7055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90543,0.0263147,-2.79303e-05,1.96928e-08,-6.1986e-12,23625.6,12.182], Tmin=(100,'K'), Tmax=(747.867,'K')), NASAPolynomial(coeffs=[4.81186,0.016118,-7.47822e-06,1.46101e-09,-1.03887e-13,23340.4,3.53847], Tmin=(747.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(F)[C]=C(7717)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {4,D} {8,S} {9,S}
4  C u1 p0 c0 {2,S} {3,D}
5  C u0 p0 c0 {2,S} {6,T}
6  C u0 p0 c0 {5,T} {10,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (290.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.4907,'amu*angstrom^2'), symmetry=1, barrier=(34.2741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4886,'amu*angstrom^2'), symmetry=1, barrier=(34.2258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69373,0.0545804,-7.68189e-05,6.28245e-08,-2.073e-11,34959.1,18.1814], Tmin=(100,'K'), Tmax=(830.634,'K')), NASAPolynomial(coeffs=[7.18289,0.0229495,-1.03126e-05,1.91362e-09,-1.30096e-13,34226.5,-6.20272], Tmin=(830.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(F)C#C(7718)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u0 p0 c0 {2,S} {6,T}
5  C u1 p0 c0 {3,D} {9,S}
6  C u0 p0 c0 {4,T} {10,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (299.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2175,525,3120,650,792.5,1650,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.50442,'amu*angstrom^2'), symmetry=1, barrier=(34.5895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50472,'amu*angstrom^2'), symmetry=1, barrier=(34.5964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0835,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72363,0.0525823,-6.50146e-05,4.54306e-08,-1.29526e-11,36072.5,18.0091], Tmin=(100,'K'), Tmax=(851.759,'K')), NASAPolynomial(coeffs=[8.4777,0.0208619,-9.14916e-06,1.70206e-09,-1.16979e-13,34922,-13.4912], Tmin=(851.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P)"""),
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
    label = '[CH]=C=C[C]=C(7719)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {3,S} {4,D} {6,S}
2 C u0 p0 c0 {3,D} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {2,D}
4 C u0 p0 c0 {1,D} {5,D}
5 C u1 p0 c0 {4,D} {9,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (587.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.7607,'amu*angstrom^2'), symmetry=1, barrier=(40.4819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0851,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63527,0.0441681,-4.52116e-05,2.40403e-08,-4.84463e-12,70808.6,17.5612], Tmin=(100,'K'), Tmax=(1402.08,'K')), NASAPolynomial(coeffs=[11.4693,0.00977197,-1.62974e-06,9.2344e-11,5.80485e-16,68674.2,-30.9832], Tmin=(1402.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=CC(F)=C[CH2](7720)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {5,S} {7,S}
4  C u0 p0 c0 {2,S} {6,D} {8,S}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u1 p0 c0 {4,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (227.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78851,0.0403144,-1.30524e-05,-1.24756e-08,7.72969e-12,27503.8,20.2707], Tmin=(100,'K'), Tmax=(1008.26,'K')), NASAPolynomial(coeffs=[11.9227,0.0181781,-7.00061e-06,1.29649e-09,-9.20879e-14,24541.8,-33.2587], Tmin=(1008.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(F)[C]=C(7721)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {6,D} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {8,S} {11,S}
5  C u1 p0 c0 {2,S} {4,D}
6  C u1 p0 c0 {2,S} {3,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (361.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72953,0.0558427,-8.14331e-05,7.45673e-08,-2.78365e-11,43602.4,20.2138], Tmin=(100,'K'), Tmax=(791.403,'K')), NASAPolynomial(coeffs=[4.70116,0.0308837,-1.52875e-05,2.97746e-09,-2.08491e-13,43443.3,8.53917], Tmin=(791.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(F)C=[CH](7722)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {6,D} {8,S}
4  C u0 p0 c0 {2,S} {5,D} {9,S}
5  C u1 p0 c0 {4,D} {10,S}
6  C u1 p0 c0 {3,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (380.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(1.0518,'amu*angstrom^2'), symmetry=1, barrier=(24.1829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05125,'amu*angstrom^2'), symmetry=1, barrier=(24.1703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89689,0.0504752,-5.24661e-05,3.18199e-08,-8.3209e-12,45824.6,19.4896], Tmin=(100,'K'), Tmax=(898.505,'K')), NASAPolynomial(coeffs=[7.19977,0.0268673,-1.30537e-05,2.57652e-09,-1.84111e-13,44871.6,-5.52612], Tmin=(898.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(F)C=C(7723)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {3,D} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {6,D}
6  C u1 p0 c0 {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (371.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.972207,'amu*angstrom^2'), symmetry=1, barrier=(22.353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970968,'amu*angstrom^2'), symmetry=1, barrier=(22.3245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15494,0.0475737,-3.7017e-05,-1.00547e-08,2.7552e-11,44698.8,19.4013], Tmin=(100,'K'), Tmax=(502.947,'K')), NASAPolynomial(coeffs=[5.84489,0.0290944,-1.43147e-05,2.81433e-09,-1.99613e-13,44190.2,2.76898], Tmin=(502.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C(F)[C]=C(7601)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,D} {4,S} {7,S}
3  C u0 p0 c0 {1,S} {2,D} {6,S}
4  C u1 p0 c0 {2,S} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (201.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,86,203,488,582,605,741,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(2.52146,'amu*angstrom^2'), symmetry=1, barrier=(57.9734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11103,'amu*angstrom^2'), symmetry=1, barrier=(57.9196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63199,0.0471134,-3.70022e-05,1.52188e-08,-2.54571e-12,24288,21.4346], Tmin=(100,'K'), Tmax=(1411.93,'K')), NASAPolynomial(coeffs=[11.3878,0.0194747,-7.63912e-06,1.35434e-09,-9.07871e-14,21533.2,-28.9968], Tmin=(1411.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cdj(Cd-F1sCd)(Cd-HH))"""),
)

species(
    label = '[CH]=CC=C([CH2])F(7724)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,D} {4,S} {7,S}
3  C u0 p0 c0 {1,S} {2,D} {5,S}
4  C u0 p0 c0 {2,S} {6,D} {8,S}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u1 p0 c0 {4,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (214.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68785,0.0420053,-1.37418e-05,-1.64592e-08,1.06533e-11,25895.6,20.1585], Tmin=(100,'K'), Tmax=(952.818,'K')), NASAPolynomial(coeffs=[13.2164,0.0149731,-4.82087e-06,8.32785e-10,-5.85966e-14,22728.8,-39.9921], Tmin=(952.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C=C[CH]F(6685)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,D} {4,S} {7,S}
3  C u0 p0 c0 {2,D} {6,S} {8,S}
4  C u1 p0 c0 {1,S} {2,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (212.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,234,589,736,816,1240,3237,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.59776,'amu*angstrom^2'), symmetry=1, barrier=(36.7357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59733,'amu*angstrom^2'), symmetry=1, barrier=(36.7258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0915,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3828,0.0512814,-4.51493e-05,2.06388e-08,-3.75853e-12,25675.4,19.9141], Tmin=(100,'K'), Tmax=(1326.52,'K')), NASAPolynomial(coeffs=[12.8001,0.0168537,-6.21922e-06,1.07376e-09,-7.12472e-14,22646.3,-38.3936], Tmin=(1326.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(C=CJC=C)"""),
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
    E0 = (88.9529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (183.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (502.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (97.2372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (167.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (336.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (324.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (237.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (174.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (236.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (232.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (504.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (321.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (271.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (194.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (203.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (264.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (281.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (202.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (245.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['C2H2(23)', 'C=C=CF(6101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['[CH]=CC([CH2])=CF(7410)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C(395)', '[CH]C=CF(270)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['C=C1C=CC1F(7714)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['C=C=C(F)C=C(7715)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[CH]=CC=C=C(6456)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(63.9786,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[CH](264)', 'C=C=CF(6101)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(323.516,'m^3/(mol*s)'), n=1.13767, Ea=(6.81221,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.011839139166081196, var=0.8671401158115553, Tref=1000.0, N=160, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', '[CH]=CC(F)=C=C(7716)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(58.9142,'m^3/(mol*s)'), n=1.72001, Ea=(9.72174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2675934934080712, var=0.0038470988293168424, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H2(23)', 'C=[C][CH]F(6102)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.81516e-05,'m^3/(mol*s)'), n=3.04336, Ea=(42.9923,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'C#CC(F)[C]=C(7717)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(22.785,'m^3/(mol*s)'), n=1.84735, Ea=(16.9899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.2364442503798854, var=2.6394824742839527, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_Sp-6R!H-4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_Sp-6R!H-4CClOS"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[CH]=CC(F)C#C(7718)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(172591,'m^3/(mol*s)'), n=0.983154, Ea=(3.98016,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[CH](264)', 'C=[C][CH]F(6102)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(408.814,'m^3/(mol*s)'), n=1.30994, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.006495679117668562, var=10.327386067777251, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_N-4R!H->O_N-4BrCClF->F',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_N-4R!H->O_N-4BrCClF->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HF(38)', '[CH]=C=C[C]=C(7719)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(297.017,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['[CH]=CC(F)=C[CH2](7720)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['C=[C]C(F)[C]=C(7721)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC(F)C=[CH](7722)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.16e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C(F)C=C(7723)'],
    products = ['[CH]=CC(F)[C]=C(7411)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['[CH2]C=C(F)[C]=C(7601)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['[CH]=CC=C([CH2])F(7724)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(113.279,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC(F)[C]=C(7411)'],
    products = ['C=[C]C=C[CH]F(6685)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(156.085,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #2917',
    isomers = [
        '[CH]=CC(F)[C]=C(7411)',
    ],
    reactants = [
        ('C2H2(23)', 'C=C=CF(6101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2917',
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

