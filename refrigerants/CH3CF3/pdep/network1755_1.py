species(
    label = 'C=[C]OC(F)(F)[C]=CF(6171)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {3,S} {8,D} {10,S}
7  C u0 p0 c0 {9,D} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {6,D}
9  C u1 p0 c0 {4,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-231.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,615,860,1140,1343,3152,2950,3100,1380,975,1025,1650,1670,1700,300,440,273.57,273.577,273.579,273.585],'cm^-1')),
        HinderedRotor(inertia=(0.279659,'amu*angstrom^2'), symmetry=1, barrier=(14.8531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279654,'amu*angstrom^2'), symmetry=1, barrier=(14.8528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279656,'amu*angstrom^2'), symmetry=1, barrier=(14.8526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.233623,0.0855846,-0.000112768,7.40549e-08,-1.9079e-11,-27746.9,29.6391], Tmin=(100,'K'), Tmax=(952.231,'K')), NASAPolynomial(coeffs=[15.809,0.0201577,-9.70418e-06,1.89859e-09,-1.34913e-13,-30713.2,-44.7406], Tmin=(952.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(C=CJO)"""),
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
    label = 'FC=C=C(F)F(1325)',
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
    label = '[C]=C(976)',
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
    label = '[O]C(F)(F)[C]=CF(547)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6 C u0 p0 c0 {3,S} {7,D} {8,S}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-300.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.813738,'amu*angstrom^2'), symmetry=1, barrier=(18.7094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97056,0.0474072,-5.58215e-05,3.38039e-08,-8.23918e-12,-36115.3,20.5538], Tmin=(100,'K'), Tmax=(990.624,'K')), NASAPolynomial(coeffs=[9.79057,0.0158305,-8.0071e-06,1.62531e-09,-1.18216e-13,-37664.6,-17.0994], Tmin=(990.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sF1sCd)) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
)

species(
    label = '[C]=CF-2(1206)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u2 p0 c0 {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (404.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([194,682,905,1196,1383,3221],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62308,0.00704726,-1.17425e-06,-1.98484e-09,8.12221e-13,48709.9,8.54799], Tmin=(100,'K'), Tmax=(1284.6,'K')), NASAPolynomial(coeffs=[5.40191,0.00467991,-2.11332e-06,4.24428e-10,-3.06823e-14,47991.2,-1.49789], Tmin=(1284.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C=[C]O[C](F)F(1702)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 C u0 p0 c0 {6,D} {7,S} {8,S}
6 C u1 p0 c0 {3,S} {5,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-119.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,2950,3100,1380,975,1025,1650,1685,370,180,180,2136.55],'cm^-1')),
        HinderedRotor(inertia=(0.270111,'amu*angstrom^2'), symmetry=1, barrier=(6.21039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270364,'amu*angstrom^2'), symmetry=1, barrier=(6.2162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16002,0.0488447,-9.07756e-05,9.20773e-08,-3.52103e-11,-14301.1,21.4717], Tmin=(100,'K'), Tmax=(833.524,'K')), NASAPolynomial(coeffs=[3.27168,0.0245568,-1.29593e-05,2.55841e-09,-1.7865e-13,-13828,20.2605], Tmin=(833.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(F1s)(F1s)(O2s-Cd)) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OC(F)(F)C1=CF(6178)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {8,D}
7  C u0 p0 c0 {4,S} {6,S} {9,D}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  C u0 p0 c0 {7,D} {11,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-615.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372531,0.0685568,-5.48918e-05,1.27983e-08,2.06542e-12,-73938,20.5052], Tmin=(100,'K'), Tmax=(1023,'K')), NASAPolynomial(coeffs=[19.9787,0.0133576,-5.42375e-06,1.06862e-09,-7.96898e-14,-79072.5,-80.0178], Tmin=(1023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCFH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'F[CH][C]=C(F)F(2842)',
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
    label = 'C=[C]OC(F)=C=CF(6294)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {7,D}
5  C u0 p0 c0 {2,S} {7,D} {11,S}
6  C u0 p0 c0 {8,D} {9,S} {10,S}
7  C u0 p0 c0 {4,D} {5,D}
8  C u1 p0 c0 {3,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (10.4442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([197,221,431,657,113,247,382,1207,3490,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,543.409,544.592,1830.38,1832.24,1832.88],'cm^-1')),
        HinderedRotor(inertia=(0.167626,'amu*angstrom^2'), symmetry=1, barrier=(35.1217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315829,'amu*angstrom^2'), symmetry=1, barrier=(6.57897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19013,0.0677269,-0.000101222,8.76934e-08,-3.09e-11,1351.63,28.327], Tmin=(100,'K'), Tmax=(777.143,'K')), NASAPolynomial(coeffs=[7.40617,0.0293228,-1.47254e-05,2.87954e-09,-2.02039e-13,579.039,1.15072], Tmin=(777.143,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.4442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCddFO) + group(Cds-CdsHH) + group(CdCddFH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=O(981)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101366,2.30728e-06,-8.97551e-09,3.68236e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35057,0.00638948,-2.69366e-06,5.42206e-10,-4.02473e-14,18240.9,-6.33612], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'C#CC(F)(F)O[C]=C(6295)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {4,S} {7,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
5  C u0 p0 c0 {7,D} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,T}
7  C u1 p0 c0 {3,S} {5,D}
8  C u0 p0 c0 {6,T} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-137.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([157,399,436,611,601,744,1127,2950,3100,1380,975,1025,1650,2175,525,1685,370,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.26873,'amu*angstrom^2'), symmetry=1, barrier=(29.1705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26782,'amu*angstrom^2'), symmetry=1, barrier=(29.1496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27021,'amu*angstrom^2'), symmetry=1, barrier=(29.2046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320858,0.0758071,-9.06589e-05,5.08731e-08,-1.04979e-11,-16383.7,25.0622], Tmin=(100,'K'), Tmax=(970.689,'K')), NASAPolynomial(coeffs=[18.5153,0.0100045,-3.14888e-06,5.06484e-10,-3.32696e-14,-20348.1,-64.4001], Tmin=(970.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO)"""),
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
    label = 'C=[C]OC(F)(F)C#CF(6296)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {8,D} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {9,T}
8  C u1 p0 c0 {4,S} {6,D}
9  C u0 p0 c0 {3,S} {7,T}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-237.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([157,399,436,611,601,744,1127,2950,3100,1380,975,1025,1650,2175,525,1685,370,239,401,1367,428.321,428.55,429.034,429.075],'cm^-1')),
        HinderedRotor(inertia=(0.316321,'amu*angstrom^2'), symmetry=1, barrier=(41.2077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109954,'amu*angstrom^2'), symmetry=1, barrier=(14.3407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109831,'amu*angstrom^2'), symmetry=1, barrier=(14.3464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.326877,0.0820025,-0.000109656,7.18983e-08,-1.83094e-11,-28388.5,27.3707], Tmin=(100,'K'), Tmax=(967.87,'K')), NASAPolynomial(coeffs=[16.3084,0.0159543,-7.294e-06,1.39171e-09,-9.75761e-14,-31482.1,-49.2088], Tmin=(967.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(CtCF) + radical(C=CJO)"""),
)

species(
    label = 'C#COC(F)(F)[C]=CF(6297)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {7,D} {10,S}
7  C u1 p0 c0 {5,S} {6,D}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-262.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,615,860,1140,1343,3152,1685,370,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13594,'amu*angstrom^2'), symmetry=1, barrier=(26.1174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13688,'amu*angstrom^2'), symmetry=1, barrier=(26.1391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13699,'amu*angstrom^2'), symmetry=1, barrier=(26.1416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00489307,0.0898313,-0.000123392,8.09573e-08,-2.05215e-11,-31395.8,24.8213], Tmin=(100,'K'), Tmax=(973.441,'K')), NASAPolynomial(coeffs=[18.3281,0.0144992,-7.31168e-06,1.45967e-09,-1.05017e-13,-34965.1,-63.1313], Tmin=(973.441,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
)

species(
    label = 'C=[C]OC(F)(F)C=[C]F(6298)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {9,D} {10,S}
7  C u0 p0 c0 {8,D} {11,S} {12,S}
8  C u1 p0 c0 {4,S} {7,D}
9  C u1 p0 c0 {3,S} {6,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-244.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193043,0.0836964,-0.000103553,6.26101e-08,-1.47505e-11,-29327.6,29.2567], Tmin=(100,'K'), Tmax=(1042.07,'K')), NASAPolynomial(coeffs=[17.4524,0.0174467,-8.19157e-06,1.60271e-09,-1.14558e-13,-32924.7,-54.7212], Tmin=(1042.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[CH]=COC(F)(F)[C]=CF(6299)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {4,S} {9,D} {10,S}
7  C u0 p0 c0 {3,S} {8,D} {11,S}
8  C u1 p0 c0 {5,S} {7,D}
9  C u1 p0 c0 {6,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-224.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,3010,987.5,1337.5,450,1655,615,860,1140,1343,3152,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.94482,'amu*angstrom^2'), symmetry=1, barrier=(21.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942483,'amu*angstrom^2'), symmetry=1, barrier=(21.6695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939782,'amu*angstrom^2'), symmetry=1, barrier=(21.6074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.28552,0.0902562,-0.000113273,6.72904e-08,-1.5303e-11,-26838.2,28.7674], Tmin=(100,'K'), Tmax=(1088.92,'K')), NASAPolynomial(coeffs=[21.1386,0.0115575,-4.86492e-06,9.19658e-10,-6.52434e-14,-31504,-76.4166], Tmin=(1088.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(F)(F)C=CF(6300)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {3,S} {6,D} {11,S}
8  C u1 p0 c0 {4,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-248.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.106517,0.0847459,-0.000104568,6.25625e-08,-1.45308e-11,-29687.2,29.4785], Tmin=(100,'K'), Tmax=(1058.65,'K')), NASAPolynomial(coeffs=[18.2254,0.0162857,-7.56737e-06,1.4779e-09,-1.05704e-13,-33523.5,-58.9675], Tmin=(1058.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=COC(F)(F)[C]=[C]F(6301)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {9,D}
9  C u1 p0 c0 {3,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-221.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.904691,'amu*angstrom^2'), symmetry=1, barrier=(20.8006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905489,'amu*angstrom^2'), symmetry=1, barrier=(20.819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905746,'amu*angstrom^2'), symmetry=1, barrier=(20.8249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.191774,0.0891276,-0.000112008,6.70468e-08,-1.54108e-11,-26478.9,28.5194], Tmin=(100,'K'), Tmax=(1075.88,'K')), NASAPolynomial(coeffs=[20.3565,0.0127309,-5.49499e-06,1.04563e-09,-7.41788e-14,-30900.4,-72.1172], Tmin=(1075.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(CdCFH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'C=[C]O[C](F)C(F)=CF(6302)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,D}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {5,D} {10,S}
8  C u0 p0 c0 {9,D} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-228.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,2950,3100,1380,975,1025,1650,1685,370,214.35,214.521,214.861,215.029],'cm^-1')),
        HinderedRotor(inertia=(0.826403,'amu*angstrom^2'), symmetry=1, barrier=(26.8019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81198,'amu*angstrom^2'), symmetry=1, barrier=(26.7905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.808223,'amu*angstrom^2'), symmetry=1, barrier=(26.8174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202136,0.0828457,-0.0001006,5.97747e-08,-1.3842e-11,-27333.2,29.5708], Tmin=(100,'K'), Tmax=(1059.84,'K')), NASAPolynomial(coeffs=[17.5246,0.017467,-8.06736e-06,1.56831e-09,-1.11788e-13,-31004.9,-55.0064], Tmin=(1059.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = 'C=C(F)OC(F)=[C][CH]F(6303)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {9,D}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  C u1 p0 c0 {3,S} {9,S} {12,S}
9  C u1 p0 c0 {6,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-227.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,293,496,537,1218,2950,3100,1380,975,1025,1650,234,589,736,816,1240,3237,1685,370,180,180,1080.62,1080.77],'cm^-1')),
        HinderedRotor(inertia=(0.233917,'amu*angstrom^2'), symmetry=1, barrier=(5.37822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233696,'amu*angstrom^2'), symmetry=1, barrier=(5.37313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23424,'amu*angstrom^2'), symmetry=1, barrier=(5.38563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410864,0.0914478,-0.000163092,1.55817e-07,-5.71185e-11,-27198.5,29.3609], Tmin=(100,'K'), Tmax=(832.616,'K')), NASAPolynomial(coeffs=[6.10273,0.0381295,-2.02431e-05,3.97315e-09,-2.76285e-13,-27246,8.35012], Tmin=(832.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCFHH) + group(Cds-CdsCsH) + group(CdCFO) + group(CdCFO) + group(Cds-CdsHH) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-F1sO2s))"""),
)

species(
    label = '[CH]=C(F)C(F)(F)O[C]=C(6304)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {3,S} {5,S} {9,D}
7  C u0 p0 c0 {8,D} {10,S} {11,S}
8  C u1 p0 c0 {4,S} {7,D}
9  C u1 p0 c0 {6,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-236.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0116391,0.0866241,-0.000108818,6.62356e-08,-1.56076e-11,-28349.4,29.3144], Tmin=(100,'K'), Tmax=(1045.71,'K')), NASAPolynomial(coeffs=[18.558,0.0156796,-7.0509e-06,1.35461e-09,-9.60296e-14,-32228.1,-60.9897], Tmin=(1045.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-236.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cdj(Cd-CsF1s)(H))"""),
)

species(
    label = '[CH]=[C]C(F)(F)OC(=C)F(6305)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {10,S} {11,S}
8  C u1 p0 c0 {5,S} {9,D}
9  C u1 p0 c0 {8,D} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-231.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,326,540,652,719,1357,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.738487,'amu*angstrom^2'), symmetry=1, barrier=(16.9793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739082,'amu*angstrom^2'), symmetry=1, barrier=(16.993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.996474,'amu*angstrom^2'), symmetry=1, barrier=(22.9109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242158,0.0888007,-0.000123382,8.77058e-08,-2.4845e-11,-27752,27.2097], Tmin=(100,'K'), Tmax=(861.919,'K')), NASAPolynomial(coeffs=[13.8461,0.0256672,-1.35102e-05,2.72352e-09,-1.9585e-13,-30097.1,-36.4001], Tmin=(861.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sO2s)(Cd-HH)) + radical(Cds_P)"""),
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
    E0 = (-138.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (392.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (378.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-130.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-94.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (207.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-105.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (97.2014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (67.5315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (50.7963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (63.6293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (65.1814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-26.1672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (106.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (180.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (16.5723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (91.0505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (12.1889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (58.9567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    products = ['CH2CO(28)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=C(976)', '[O]C(F)(F)[C]=CF(547)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=CF-2(1206)', 'C=[C]O[C](F)F(1702)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    products = ['C=C1OC(F)(F)C1=CF(6178)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CO(28)', 'F[CH][C]=C(F)F(2842)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000537,'m^3/(mol*s)'), n=2.41, Ea=(74.4905,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Sp-5R!H=4R!H_Ext-2R!H-R_Sp-2R!H=1R!H_N-Sp-8R!H-2R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Sp-5R!H=4R!H_Ext-2R!H-R_Sp-2R!H=1R!H_N-Sp-8R!H-2R!H"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'C=[C]OC(F)=C=CF(6294)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(31.4053,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][C]=O(981)', 'FC=C=C(F)F(1325)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(5.69712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'C#CC(F)(F)O[C]=C(6295)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(68.8176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C=[C]OC(F)(F)C#CF(6296)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'C#COC(F)(F)[C]=CF(6297)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(125.384,'m^3/(mol*s)'), n=2.03877, Ea=(8.36527,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_N-1COS->O_N-Sp-2CS=1CCSS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=O(981)', 'F[CH][C]=C(F)F(2842)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(11.2564,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    products = ['C=[C]OC(F)(F)C=[C]F(6298)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=COC(F)(F)[C]=CF(6299)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    products = ['[CH]=[C]OC(F)(F)C=CF(6300)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=COC(F)(F)[C]=[C]F(6301)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_single;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_single;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]O[C](F)C(F)=CF(6302)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(152.123,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C(F)OC(F)=[C][CH]F(6303)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(225.312,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    products = ['[CH]=C(F)C(F)(F)O[C]=C(6304)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(151.145,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]C(F)(F)OC(=C)F(6305)'],
    products = ['C=[C]OC(F)(F)[C]=CF(6171)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(197.927,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1755',
    isomers = [
        'C=[C]OC(F)(F)[C]=CF(6171)',
    ],
    reactants = [
        ('CH2CO(28)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1755',
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

