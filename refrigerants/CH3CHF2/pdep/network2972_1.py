species(
    label = 'C=[C]C(F)O[CH]F(7462)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
5  C u1 p0 c0 {2,S} {3,S} {11,S}
6  C u0 p0 c0 {7,D} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-131.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,580,1155,1237,1373,3147,2950,3100,1380,975,1025,1650,1685,370,283.204,283.214,283.236],'cm^-1')),
        HinderedRotor(inertia=(0.166155,'amu*angstrom^2'), symmetry=1, barrier=(9.45681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166258,'amu*angstrom^2'), symmetry=1, barrier=(9.45633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531959,'amu*angstrom^2'), symmetry=1, barrier=(30.2762,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3439,0.0636573,-8.38971e-05,6.26082e-08,-1.93489e-11,-15769.4,23.9649], Tmin=(100,'K'), Tmax=(781.692,'K')), NASAPolynomial(coeffs=[8.49612,0.0270587,-1.36679e-05,2.71327e-09,-1.93462e-13,-16887.6,-8.77888], Tmin=(781.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsFHHO) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(F1s)(O2s-Cs)(H)) + radical(Cdj(Cs-F1sO2sH)(Cd-HH))"""),
)

species(
    label = 'CHFO(47)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.294,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]F(230)',
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
    label = 'C=[C]C([O])F(714)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (96.8427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.861269,'amu*angstrom^2'), symmetry=1, barrier=(19.8023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4047,0.035919,-3.5484e-05,1.85546e-08,-3.95241e-12,11704.2,16.1894], Tmin=(100,'K'), Tmax=(1119.75,'K')), NASAPolynomial(coeffs=[8.47789,0.0142241,-6.42174e-06,1.25173e-09,-8.9283e-14,10344.2,-13.7971], Tmin=(1119.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.8427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCdH)) + radical(Cdj(Cs-F1sO2sH)(Cd-HH))"""),
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
    label = 'F[CH]O[CH]F(639)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u1 p0 c0 {1,S} {3,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-238.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([505,655,1117,1193,1185,1289,1311,1435,3104,3190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.511832,'amu*angstrom^2'), symmetry=1, barrier=(11.768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59482,'amu*angstrom^2'), symmetry=1, barrier=(13.6761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02831,0.0488021,-8.68646e-05,7.88728e-08,-2.74521e-11,-28640,16.0059], Tmin=(100,'K'), Tmax=(838.02,'K')), NASAPolynomial(coeffs=[7.11357,0.0145015,-7.5195e-06,1.47268e-09,-1.02147e-13,-29140.2,-5.52802], Tmin=(838.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(F1s)(O2s-Cs)(H)) + radical(Csj(F1s)(O2s-Cs)(H))"""),
)

species(
    label = 'C=C1C(F)OC1F(7465)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-455.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83896,0.0348232,1.1011e-05,-3.82546e-08,1.61072e-11,-54664.1,14.8939], Tmin=(100,'K'), Tmax=(1037.92,'K')), NASAPolynomial(coeffs=[13.6288,0.0189044,-8.64162e-06,1.76833e-09,-1.32661e-13,-58701.4,-50.0833], Tmin=(1037.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-455.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsCFHO) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=C=C(F)OCF(7857)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u0 p0 c0 {7,D} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-387.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15312,0.057581,-5.48523e-05,2.57258e-08,-4.75127e-12,-46517.8,21.2443], Tmin=(100,'K'), Tmax=(1310.23,'K')), NASAPolynomial(coeffs=[14.8841,0.0156614,-6.86084e-06,1.30678e-09,-9.19351e-14,-50115.9,-48.7097], Tmin=(1310.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-387.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFHHO) + group(CdCddFO) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C=CO[CH]F(6553)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {3,S} {4,S}
3  C u0 p0 c0 {2,S} {6,D} {7,S}
4  C u1 p0 c0 {1,S} {2,S} {8,S}
5  C u0 p0 c0 {6,D} {9,S} {10,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-24.0321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,580,1155,1237,1373,3147,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21561,'amu*angstrom^2'), symmetry=1, barrier=(27.9493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2126,'amu*angstrom^2'), symmetry=1, barrier=(27.88,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0722,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29471,0.0486343,-2.64137e-05,-1.30246e-08,1.19695e-11,-2782.98,19.7746], Tmin=(100,'K'), Tmax=(933.355,'K')), NASAPolynomial(coeffs=[17.9506,0.00483432,-3.47903e-07,1.77382e-11,-4.16392e-15,-7093.51,-65.8676], Tmin=(933.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.0321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFHHO) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Csj(F1s)(O2s-Cd)(H))"""),
)

species(
    label = '[O][CH]F(516)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53591,0.00799828,-5.22385e-07,-4.56598e-09,2.08746e-12,-2348.33,9.31393], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97189,0.00388217,-1.62977e-06,3.36434e-10,-2.54185e-14,-3160.19,-3.94933], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
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
    label = 'C=C=C(F)O[CH]F(7858)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {7,D}
5  C u1 p0 c0 {2,S} {3,S} {8,S}
6  C u0 p0 c0 {7,D} {9,S} {10,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-189.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([197,221,431,657,580,1155,1237,1373,3147,2950,3100,1380,975,1025,1650,540,610,2055,408.681,408.951,410.396,412.552],'cm^-1')),
        HinderedRotor(inertia=(0.159183,'amu*angstrom^2'), symmetry=1, barrier=(18.8294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(18.8025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25406,0.0596299,-6.7009e-05,3.72598e-08,-8.15132e-12,-22748.2,21.6289], Tmin=(100,'K'), Tmax=(1113.78,'K')), NASAPolynomial(coeffs=[13.6444,0.0151309,-7.07879e-06,1.38742e-09,-9.93012e-14,-25508.2,-39.4826], Tmin=(1113.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFHHO) + group(CdCddFO) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Csj(F1s)(O2s-Cd)(H))"""),
)

species(
    label = 'C#CC(F)O[CH]F(7859)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
5  C u1 p0 c0 {2,S} {3,S} {9,S}
6  C u0 p0 c0 {4,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-208.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([242,464,653,889,1056,1342,3790,580,1155,1237,1373,3147,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.45316,'amu*angstrom^2'), symmetry=1, barrier=(56.4029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409029,'amu*angstrom^2'), symmetry=1, barrier=(12.456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.45467,'amu*angstrom^2'), symmetry=1, barrier=(56.4377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15596,0.0672717,-0.000102456,8.42827e-08,-2.75013e-11,-24963,21.0623], Tmin=(100,'K'), Tmax=(822.983,'K')), NASAPolynomial(coeffs=[9.37479,0.0220576,-1.04464e-05,1.97192e-09,-1.34969e-13,-26137.4,-15.9038], Tmin=(822.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsFHHO) + group(Ct-CtCs) + group(Ct-CtH) + radical(Csj(F1s)(O2s-Cs)(H))"""),
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
    label = '[CH2]C=C(F)O[CH]F(7860)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {7,S}
4  C u0 p0 c0 {5,D} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {4,D}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {3,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-234.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03236,0.0604461,-6.05761e-05,2.99741e-08,-5.81951e-12,-28074,23.3552], Tmin=(100,'K'), Tmax=(1253.26,'K')), NASAPolynomial(coeffs=[15.1748,0.0153075,-6.5503e-06,1.23512e-09,-8.66214e-14,-31618.9,-48.0666], Tmin=(1253.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-234.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(CsFHHO) + group(Cds-CdsCsH) + group(CdCFO) + radical(Allyl_P) + radical(Csj(F1s)(O2s-Cd)(H))"""),
)

species(
    label = '[CH]=CC(F)O[CH]F(7861)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {4,S} {6,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u1 p0 c0 {2,S} {3,S} {10,S}
7  C u1 p0 c0 {5,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-127.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([236,527,855,1015,1182,1348,3236,3010,987.5,1337.5,450,1655,580,1155,1237,1373,3147,3120,650,792.5,1650,224.668,224.688,224.707],'cm^-1')),
        HinderedRotor(inertia=(0.245723,'amu*angstrom^2'), symmetry=1, barrier=(8.80371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245716,'amu*angstrom^2'), symmetry=1, barrier=(8.80365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931802,'amu*angstrom^2'), symmetry=1, barrier=(33.3759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23583,0.0663319,-9.43728e-05,7.72223e-08,-2.61271e-11,-15207,23.5671], Tmin=(100,'K'), Tmax=(717.227,'K')), NASAPolynomial(coeffs=[8.09473,0.0280778,-1.43646e-05,2.85062e-09,-2.02517e-13,-16190.8,-7.24316], Tmin=(717.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsFHHO) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Csj(F1s)(O2s-Cs)(H)) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=C(F)OCF(7862)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {7,S} {10,S} {11,S}
7  C u1 p0 c0 {5,D} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-182.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25182,0.0585138,-5.76613e-05,2.85782e-08,-5.65014e-12,-21823.9,22.6355], Tmin=(100,'K'), Tmax=(1217.06,'K')), NASAPolynomial(coeffs=[13.3947,0.018605,-8.47473e-06,1.63546e-09,-1.15779e-13,-24779.7,-38.3324], Tmin=(1217.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(CsFHHO) + group(Cds-CdsCsH) + group(CdCFO) + radical(Allyl_P) + radical(Cdj(Cs-HHH)(Cd-F1sO2s))"""),
)

species(
    label = '[CH]=[C]C(F)OCF(7863)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
5  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
6  C u1 p0 c0 {4,S} {7,D}
7  C u1 p0 c0 {6,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-83.4642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,548,1085,1183,1302,1466,1520,3060,3119,1685,370,3120,650,792.5,1650,243.075,243.11],'cm^-1')),
        HinderedRotor(inertia=(0.724782,'amu*angstrom^2'), symmetry=1, barrier=(30.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364131,'amu*angstrom^2'), symmetry=1, barrier=(15.2697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244904,'amu*angstrom^2'), symmetry=1, barrier=(10.2662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47285,0.0602347,-7.04582e-05,4.52811e-08,-1.20724e-11,-9951.39,23.6524], Tmin=(100,'K'), Tmax=(898.063,'K')), NASAPolynomial(coeffs=[9.20466,0.0257965,-1.29366e-05,2.57998e-09,-1.85209e-13,-11340.1,-12.8177], Tmin=(898.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.4642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsFHHO) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sO2sH)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=COC(F)F(6675)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {4,S} {5,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
5  C u0 p0 c0 {3,S} {7,D} {9,S}
6  C u1 p0 c0 {7,S} {10,S} {11,S}
7  C u1 p0 c0 {5,D} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-241.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190791,0.0665036,-6.8669e-05,3.32103e-08,-5.98732e-12,-28884.6,25.6576], Tmin=(100,'K'), Tmax=(1548.01,'K')), NASAPolynomial(coeffs=[20.2613,0.00538431,-4.74676e-07,-2.15642e-11,3.46456e-15,-33989.2,-76.3584], Tmin=(1548.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(CsFFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C(F)[CH]O[CH]F(7864)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {7,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u1 p0 c0 {3,S} {4,S} {8,S}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {3,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-239.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45408,0.0633706,-9.87445e-05,9.22202e-08,-3.43342e-11,-28717.9,24.2865], Tmin=(100,'K'), Tmax=(807.488,'K')), NASAPolynomial(coeffs=[4.79854,0.0329871,-1.66384e-05,3.24344e-09,-2.26466e-13,-28807.6,11.6557], Tmin=(807.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(CsFHHO) + group(CdCsCdF) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Csj(F1s)(O2s-Cs)(H))"""),
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
    E0 = (-101.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (342.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (391.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-93.2432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-23.2804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-67.5967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (143.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (35.0422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (52.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (38.2888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (206.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (265.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (81.0539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (8.55713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (63.3819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-20.0813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (74.8374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (37.2521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['CHFO(47)', 'C=C=CF(6101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]F(230)', 'C=[C]C([O])F(714)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C(395)', 'F[CH]O[CH]F(639)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['C=C1C(F)OC1F(7465)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_1H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['C=C=C(F)OCF(7857)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHFO(47)', 'C=[C][CH]F(6102)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(100.231,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'C=C=CO[CH]F(6553)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(64.6867,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][CH]F(516)', 'C=C=CF(6101)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(14.9195,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C=C=C(F)O[CH]F(7858)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'C#CC(F)O[CH]F(7859)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(172591,'m^3/(mol*s)'), n=0.983154, Ea=(4.50987,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_4R!H->C_Ext-4C-R_N-5R!H-inRing_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]F(516)', 'C=[C][CH]F(6102)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHF(40)', 'C=[C]C([O])F(714)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['[CH2]C=C(F)O[CH]F(7860)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC(F)O[CH]F(7861)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['[CH2][C]=C(F)OCF(7862)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3466.35,'s^-1'), n=2.68499, Ea=(164.909,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_1H;XH_out] for rate rule [R3H_SS_O;C_rad_out_1H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]C(F)OCF(7863)'],
    products = ['C=[C]C(F)O[CH]F(7462)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_1H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['[CH2][C]=COC(F)F(6675)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(176.365,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]C(F)O[CH]F(7462)'],
    products = ['C=C(F)[CH]O[CH]F(7864)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(138.78,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #2972',
    isomers = [
        'C=[C]C(F)O[CH]F(7462)',
    ],
    reactants = [
        ('CHFO(47)', 'C=C=CF(6101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2972',
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

