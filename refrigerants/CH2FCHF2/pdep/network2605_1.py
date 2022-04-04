species(
    label = 'C=[C]OC(F)[C]=C(F)F(9711)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u0 p0 c0 {9,D} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {6,D}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-195.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,562,600,623,1070,1265,2950,3100,1380,975,1025,1650,1670,1700,300,440,193.086,193.27,193.383,194.632],'cm^-1')),
        HinderedRotor(inertia=(0.816749,'amu*angstrom^2'), symmetry=1, barrier=(21.7076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.820046,'amu*angstrom^2'), symmetry=1, barrier=(21.6997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.811631,'amu*angstrom^2'), symmetry=1, barrier=(21.7079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.034916,0.0813342,-9.34304e-05,5.13994e-08,-1.08849e-11,-23401.1,31.0502], Tmin=(100,'K'), Tmax=(1164.9,'K')), NASAPolynomial(coeffs=[19.9452,0.0129657,-5.39342e-06,1.01549e-09,-7.18268e-14,-28039.8,-68.0444], Tmin=(1164.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(C=CJO)"""),
)

species(
    label = 'CH2CO(28)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'FC=C=C(F)F(5101)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {6,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-364.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([113,247,382,1207,3490,94,120,354,641,825,1294,540,610,2055,2775.82],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2930.43,'J/mol'), sigma=(4.61545,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=457.73 K, Pc=67.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78631,0.0191407,1.84311e-05,-5.99612e-08,3.72113e-11,-43861.2,10.8341], Tmin=(10,'K'), Tmax=(620.754,'K')), NASAPolynomial(coeffs=[5.19832,0.0213405,-1.41861e-05,4.38954e-09,-5.1374e-13,-44254.2,2.94178], Tmin=(620.754,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-364.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FCDCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[C]=C(4368)',
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
    label = '[O]C(F)[C]=C(F)F(9673)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-270.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,562,600,623,1070,1265,1685,370,180,180,758.456],'cm^-1')),
        HinderedRotor(inertia=(0.277283,'amu*angstrom^2'), symmetry=1, barrier=(6.37527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06534,0.0458699,-5.49839e-05,3.47006e-08,-8.89665e-12,-32523.4,21.7549], Tmin=(100,'K'), Tmax=(939.968,'K')), NASAPolynomial(coeffs=[8.9134,0.0167286,-8.48099e-06,1.71921e-09,-1.24821e-13,-33810.8,-10.8591], Tmin=(939.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-270.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCdH)) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
)

species(
    label = '[C]=C(F)F(3536)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u2 p0 c0 {3,D}
"""),
    E0 = (195.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([182,240,577,636,1210,1413],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28657,0.0148134,-1.42116e-05,6.41427e-09,-1.12277e-12,23593,10.1566], Tmin=(100,'K'), Tmax=(1382.13,'K')), NASAPolynomial(coeffs=[7.28631,0.0032378,-1.64877e-06,3.54597e-10,-2.66861e-14,22487.4,-10.4342], Tmin=(1382.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C=[C]O[CH]F(6301)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {5,S}
3 C u1 p0 c0 {1,S} {2,S} {8,S}
4 C u0 p0 c0 {5,D} {6,S} {7,S}
5 C u1 p0 c0 {2,S} {4,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (105.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,2950,3100,1380,975,1025,1650,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.210862,'amu*angstrom^2'), symmetry=1, barrier=(4.84813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210859,'amu*angstrom^2'), symmetry=1, barrier=(4.84807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27824,0.0459986,-8.51427e-05,8.6634e-08,-3.30932e-11,12716,19.0821], Tmin=(100,'K'), Tmax=(840.592,'K')), NASAPolynomial(coeffs=[2.98224,0.02393,-1.23595e-05,2.41886e-09,-1.68005e-13,13259,19.7417], Tmin=(840.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(F1s)(O2s-Cd)(H)) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OC(F)C1=C(F)F(9717)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u0 p0 c0 {4,S} {6,S} {9,D}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-582.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.527557,0.0642912,-4.34623e-05,2.08299e-09,5.5195e-12,-69892.6,19.3622], Tmin=(100,'K'), Tmax=(1017.29,'K')), NASAPolynomial(coeffs=[19.5126,0.0142064,-5.83304e-06,1.16008e-09,-8.70587e-14,-75026.4,-78.8027], Tmin=(1017.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFF) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'C=COC(F)=C=C(F)F(10829)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,D} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {9,D}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-424.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86103,0.0716542,-8.0508e-05,4.68511e-08,-1.09963e-11,-50910.5,27.5347], Tmin=(100,'K'), Tmax=(1027.37,'K')), NASAPolynomial(coeffs=[12.9237,0.0246887,-1.19363e-05,2.35438e-09,-1.6844e-13,-53389,-30.9866], Tmin=(1027.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-424.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCddFO) + group(Cds-CdsHH) + group(CdCddFF) + group(Cdd-CdsCds)"""),
)

species(
    label = 'F[CH][C]=C(F)F(5078)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u1 p0 c0 {4,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-200.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,589,736,816,1240,3237,562,600,623,1070,1265,1685,370,715.793],'cm^-1')),
        HinderedRotor(inertia=(0.35565,'amu*angstrom^2'), symmetry=1, barrier=(8.17708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37745,0.0377452,-4.40203e-05,2.68135e-08,-6.58286e-12,-24077.8,17.1544], Tmin=(100,'K'), Tmax=(983.929,'K')), NASAPolynomial(coeffs=[8.46193,0.0130101,-6.31222e-06,1.26456e-09,-9.14162e-14,-25275.2,-12.1012], Tmin=(983.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
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
    label = 'C=[C]OC=C=C(F)F(10830)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {7,D}
6  C u0 p0 c0 {8,D} {10,S} {11,S}
7  C u0 p0 c0 {4,D} {5,D}
8  C u1 p0 c0 {3,S} {6,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-18.5363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,443.753,444.644,445.533],'cm^-1')),
        HinderedRotor(inertia=(0.0882888,'amu*angstrom^2'), symmetry=1, barrier=(12.3514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149064,'amu*angstrom^2'), symmetry=1, barrier=(20.7545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06684,0.062999,-6.90395e-05,3.84512e-08,-8.46075e-12,-2122.44,28.2117], Tmin=(100,'K'), Tmax=(1107.52,'K')), NASAPolynomial(coeffs=[13.6489,0.0175566,-7.49312e-06,1.40344e-09,-9.7939e-14,-4909.4,-33.7741], Tmin=(1107.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.5363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=O(4373)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
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
    label = 'C=[C]OC(F)=C=C(F)F(10831)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {8,D}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u0 p0 c0 {9,D} {10,S} {11,S}
8  C u0 p0 c0 {5,D} {6,D}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-184.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([197,221,431,657,94,120,354,641,825,1294,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,180,1137.9,1138.12,1138.23],'cm^-1')),
        HinderedRotor(inertia=(0.129351,'amu*angstrom^2'), symmetry=1, barrier=(2.97403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129191,'amu*angstrom^2'), symmetry=1, barrier=(2.97036,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817862,0.0765683,-0.000119414,1.03088e-07,-3.57292e-11,-22078.5,30.8077], Tmin=(100,'K'), Tmax=(782.11,'K')), NASAPolynomial(coeffs=[8.8763,0.0286615,-1.46978e-05,2.88653e-09,-2.02584e-13,-23134.4,-4.78055], Tmin=(782.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCddFO) + group(Cds-CdsHH) + group(CdCddFF) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC(F)C#CF(10832)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {4,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {7,D} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,T}
7  C u1 p0 c0 {3,S} {5,D}
8  C u0 p0 c0 {2,S} {6,T}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-2.42014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([242,464,653,889,1056,1342,3790,2950,3100,1380,975,1025,1650,2175,525,1685,370,239,401,1367,396.797,396.797,396.798,396.798],'cm^-1')),
        HinderedRotor(inertia=(0.189487,'amu*angstrom^2'), symmetry=1, barrier=(21.1711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189487,'amu*angstrom^2'), symmetry=1, barrier=(21.1711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.678292,'amu*angstrom^2'), symmetry=1, barrier=(75.785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72516,0.0697238,-8.13155e-05,4.69736e-08,-1.0572e-11,-171.111,25.7286], Tmin=(100,'K'), Tmax=(1091.45,'K')), NASAPolynomial(coeffs=[15.6116,0.015167,-6.33676e-06,1.17578e-09,-8.1806e-14,-3420.66,-47.3925], Tmin=(1091.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.42014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(CtCF) + radical(C=CJO)"""),
)

species(
    label = 'C#COC(F)[C]=C(F)F(10833)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u1 p0 c0 {5,S} {6,D}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {11,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-226.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,562,600,623,1070,1265,1685,370,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35008,'amu*angstrom^2'), symmetry=1, barrier=(31.0411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35219,'amu*angstrom^2'), symmetry=1, barrier=(31.0895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34565,'amu*angstrom^2'), symmetry=1, barrier=(30.9391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.230574,0.0858782,-0.000105007,5.9439e-08,-1.278e-11,-27048.8,26.3307], Tmin=(100,'K'), Tmax=(1155.18,'K')), NASAPolynomial(coeffs=[22.4514,0.00733793,-3.02246e-06,5.82349e-10,-4.24546e-14,-32289.2,-86.3686], Tmin=(1155.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFF) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
)

species(
    label = 'C=[C]O[C](F)C=C(F)F(10834)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  C u0 p0 c0 {6,S} {7,D} {10,S}
6  C u1 p0 c0 {1,S} {4,S} {5,S}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  C u0 p0 c0 {9,D} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-256.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335141,0.0781567,-8.86426e-05,4.89868e-08,-1.05642e-11,-30677.3,29.7403], Tmin=(100,'K'), Tmax=(1134.57,'K')), NASAPolynomial(coeffs=[17.657,0.0170877,-7.9043e-06,1.54558e-09,-1.10693e-13,-34607.8,-56.0145], Tmin=(1134.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC(F)[C]=C(F)F(10835)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {4,S} {9,D} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {8,D}
8  C u1 p0 c0 {5,S} {7,D}
9  C u1 p0 c0 {6,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-188.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,3010,987.5,1337.5,450,1655,562,600,623,1070,1265,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17966,'amu*angstrom^2'), symmetry=1, barrier=(27.1228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17958,'amu*angstrom^2'), symmetry=1, barrier=(27.1208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17759,'amu*angstrom^2'), symmetry=1, barrier=(27.0751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.708646,0.0884717,-0.000101781,5.38436e-08,-1.06754e-11,-22482.3,30.9964], Tmin=(100,'K'), Tmax=(1341.76,'K')), NASAPolynomial(coeffs=[25.5645,0.00393763,-3.35822e-07,-9.90382e-12,1.36924e-15,-28973.7,-101.397], Tmin=(1341.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s)) + radical(Cds_P)"""),
)

species(
    label = 'C=CO[C](F)[C]=C(F)F(10836)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,D} {10,S}
6  C u1 p0 c0 {1,S} {4,S} {9,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-230.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.635324,0.0860942,-9.68698e-05,5.01971e-08,-9.74904e-12,-27517.7,31.0202], Tmin=(100,'K'), Tmax=(1376.8,'K')), NASAPolynomial(coeffs=[25.2305,0.00429269,-4.99174e-07,2.27481e-11,-9.81673e-16,-34009.4,-99.7479], Tmin=(1376.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFF) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sF1s))"""),
)

species(
    label = '[CH]=[C]OC(F)C=C(F)F(10837)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {6,D}
8  C u1 p0 c0 {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-214.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23085,0.0808615,-9.45402e-05,5.37188e-08,-1.18824e-11,-25640.5,29.8303], Tmin=(100,'K'), Tmax=(1108.52,'K')), NASAPolynomial(coeffs=[18.0103,0.0167065,-7.72906e-06,1.51072e-09,-1.08199e-13,-29582.3,-57.7769], Tmin=(1108.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFF) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]O[CH]C(F)=C(F)F(10838)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,D}
6  C u1 p0 c0 {4,S} {5,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,D}
8  C u0 p0 c0 {9,D} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-287.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.588121,0.0789037,-9.98356e-05,6.5662e-08,-1.72753e-11,-34423.8,29.6882], Tmin=(100,'K'), Tmax=(925.265,'K')), NASAPolynomial(coeffs=[13.1283,0.0246917,-1.19499e-05,2.33946e-09,-1.66125e-13,-36744.5,-29.837], Tmin=(925.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(CdCsCdF) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(C=CJO)"""),
)

species(
    label = 'C=C(F)O[CH][C]=C(F)F(10839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u1 p0 c0 {4,S} {9,S} {12,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,D}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-284.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11093,0.0669686,-6.81351e-05,3.61471e-08,-7.89416e-12,-34094,28.6519], Tmin=(100,'K'), Tmax=(1084.7,'K')), NASAPolynomial(coeffs=[11.765,0.0276804,-1.3805e-05,2.75573e-09,-1.98241e-13,-36405.3,-23.6143], Tmin=(1084.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(CdCFO) + group(CdCFF) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cdj(Cs-O2sHH)(Cd-F1sF1s))"""),
)

species(
    label = 'C=[C]OC(F)C(F)=[C]F(10840)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {9,D}
7  C u0 p0 c0 {8,D} {11,S} {12,S}
8  C u1 p0 c0 {4,S} {7,D}
9  C u1 p0 c0 {3,S} {6,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-173.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([236,527,855,1015,1182,1348,3236,246,474,533,1155,2950,3100,1380,975,1025,1650,1685,370,167,640,1190,402.361,402.439,402.471,402.482,402.484],'cm^-1')),
        HinderedRotor(inertia=(0.119127,'amu*angstrom^2'), symmetry=1, barrier=(13.693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245843,'amu*angstrom^2'), symmetry=1, barrier=(28.2567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119117,'amu*angstrom^2'), symmetry=1, barrier=(13.6953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.208941,0.0836535,-0.000103658,6.3193e-08,-1.50422e-11,-20674.6,29.7063], Tmin=(100,'K'), Tmax=(1030.86,'K')), NASAPolynomial(coeffs=[17.0631,0.0182552,-8.49723e-06,1.65201e-09,-1.17561e-13,-24149.4,-52.1175], Tmin=(1030.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cdj(Cd-CsF1s)(F1s))"""),
)

species(
    label = 'C=C(F)OC(F)[C]=[C]F(10841)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-172.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([162,485,641,721,926,1292,1380,3063,326,540,652,719,1357,2950,3100,1380,975,1025,1650,1685,370,167,640,1190,335.671,335.676,335.684],'cm^-1')),
        HinderedRotor(inertia=(0.319941,'amu*angstrom^2'), symmetry=1, barrier=(25.5824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319909,'amu*angstrom^2'), symmetry=1, barrier=(25.582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129627,'amu*angstrom^2'), symmetry=1, barrier=(10.3644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576715,0.0784556,-9.23556e-05,5.47664e-08,-1.29717e-11,-20624.1,27.9224], Tmin=(100,'K'), Tmax=(1022.36,'K')), NASAPolynomial(coeffs=[14.5548,0.0237659,-1.21149e-05,2.44233e-09,-1.76783e-13,-23482.2,-39.823], Tmin=(1022.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFO) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sO2sH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    E0 = (-112.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (412.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (384.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-104.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-19.1911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-78.7333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (183.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-101.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (110.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (195.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (77.4672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (52.0887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (70.2963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (0.503713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (29.6919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (133.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (33.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (72.9293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (81.2481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (97.7074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['CH2CO(28)', 'FC=C=C(F)F(5101)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=C(4368)', '[O]C(F)[C]=C(F)F(9673)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C(F)F(3536)', 'C=[C]O[CH]F(6301)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['C=C1OC(F)C1=C(F)F(9717)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['C=COC(F)=C=C(F)F(10829)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2CO(28)', 'F[CH][C]=C(F)F(5078)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(99.2387,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'C=[C]OC=C=C(F)F(10830)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(45.5323,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=O(4373)', 'FC=C=C(F)F(5101)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(19.8755,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C=[C]OC(F)=C=C(F)F(10831)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.91905e-05,'m^3/(mol*s)'), n=3.56163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2081933962573252, var=1.209330187488209, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H_Ext-1COS-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'C=[C]OC(F)C#CF(10832)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(41.4438,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'C#COC(F)[C]=C(F)F(10833)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(125.384,'m^3/(mol*s)'), n=2.03877, Ea=(8.36527,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(4373)', 'F[CH][C]=C(F)F(5078)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(9.05742,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['C=[C]O[C](F)C=C(F)F(10834)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=COC(F)[C]=C(F)F(10835)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['C=CO[C](F)[C]=C(F)F(10836)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['[CH]=[C]OC(F)C=C(F)F(10837)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['C=[C]O[CH]C(F)=C(F)F(10838)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(145.84,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    products = ['C=C(F)O[CH][C]=C(F)F(10839)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(185.214,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]OC(F)C(F)=[C]F(10840)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(170.768,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(F)OC(F)[C]=[C]F(10841)'],
    products = ['C=[C]OC(F)[C]=C(F)F(9711)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(186.678,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #2605',
    isomers = [
        'C=[C]OC(F)[C]=C(F)F(9711)',
    ],
    reactants = [
        ('CH2CO(28)', 'FC=C=C(F)F(5101)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2605',
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

