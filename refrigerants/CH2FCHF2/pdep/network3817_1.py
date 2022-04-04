species(
    label = '[O]C(F)(F)C1(F)C=C[C]1F(13756)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {6,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-503.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,346,659,817,1284,2750,3150,900,1100,180,180,180,449.171,1210.52,1212.61,1213.04],'cm^-1')),
        HinderedRotor(inertia=(1.4844,'amu*angstrom^2'), symmetry=1, barrier=(34.1294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705845,0.0741183,-7.87559e-05,4.21224e-08,-9.05244e-12,-60454.5,27.3914], Tmin=(100,'K'), Tmax=(1117.2,'K')), NASAPolynomial(coeffs=[14.532,0.0246161,-1.2293e-05,2.46253e-09,-1.77714e-13,-63543.8,-40.8442], Tmin=(1117.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-503.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sF1s)) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'CF2O(49)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-618.61,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CC=C1F(304)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {3,S} {6,D} {8,S}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (40.9583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,247.893,1569.59,1569.59,2285.13],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89822,0.00612407,0.000102608,-2.16319e-07,1.34206e-10,4932.65,9.22518], Tmin=(10,'K'), Tmax=(550.566,'K')), NASAPolynomial(coeffs=[3.92734,0.0284965,-1.98739e-05,6.49561e-09,-8.00611e-13,4587.15,5.99356], Tmin=(550.566,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(40.9583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCCDC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=C(F)[CH][CH]1(923)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {8,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (89.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,206,336,431,607,515,611,528,696,1312,1446,180,692.409,1476.48,1721.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3235.35,'J/mol'), sigma=(5.2451,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=505.35 K, Pc=50.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19216,0.0393646,-3.74e-05,1.78194e-08,-3.40083e-12,10888.1,15.3831], Tmin=(100,'K'), Tmax=(1250.71,'K')), NASAPolynomial(coeffs=[10.2921,0.0134598,-6.33194e-06,1.25922e-09,-9.07053e-14,8861.91,-25.5064], Tmin=(1250.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(Cd-CdF1s)(H)_ring) + radical(Csj(Cs-CdHH)(Cd-CdF1s)(H)_ring) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
)

species(
    label = '[O]C(F)(F)[C](F)C1(F)C=C1(13893)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
8  C u1 p0 c0 {4,S} {6,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {6,S} {9,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-345.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,253,525,597,667,842,1178,1324,212,367,445,1450,2750,3150,900,1100,180,180,529.869,532.166,537.604,537.966,542.73],'cm^-1')),
        HinderedRotor(inertia=(0.0109077,'amu*angstrom^2'), symmetry=1, barrier=(2.14077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0110989,'amu*angstrom^2'), symmetry=1, barrier=(2.25038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.391778,0.08621,-0.000123344,9.49935e-08,-2.96115e-11,-41382.7,30.144], Tmin=(100,'K'), Tmax=(781.548,'K')), NASAPolynomial(coeffs=[11.3544,0.0301001,-1.56485e-05,3.12328e-09,-2.22693e-13,-43096.1,-20.0415], Tmin=(781.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-345.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs(F)(C)-Cd) + radical(O2sj(Cs-CsF1sF1s)) + radical(CsCsCsF1s)"""),
)

species(
    label = '[O]C(F)(F)C1=CC(F)[C]1F(13894)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u1 p0 c0 {4,S} {6,S} {8,S}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-512.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476071,0.070311,-6.47563e-05,2.85169e-08,-4.92713e-12,-61486.3,28.8456], Tmin=(100,'K'), Tmax=(1395.18,'K')), NASAPolynomial(coeffs=[18.6883,0.0180963,-8.61873e-06,1.69232e-09,-1.20487e-13,-66568.2,-65.0828], Tmin=(1395.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + group(CsCCFH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-F1sF1sCd)) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)(F)C1C=C(F)[C]1F(13758)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {6,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u0 p0 c0 {4,S} {8,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-511.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,351,323,533,609,664,892,1120,1201,346,659,817,1284,271,519,563,612,1379,180,180,180,1049.28,1049.44,1049.46,1049.46,1049.52],'cm^-1')),
        HinderedRotor(inertia=(0.00832428,'amu*angstrom^2'), symmetry=1, barrier=(6.5051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3854.34,'J/mol'), sigma=(6.18801,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.04 K, Pc=36.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292476,0.071027,-6.70042e-05,3.02261e-08,-5.29299e-12,-61359.7,30.3135], Tmin=(100,'K'), Tmax=(1392.74,'K')), NASAPolynomial(coeffs=[19.8867,0.0147518,-6.39508e-06,1.21421e-09,-8.5296e-14,-66817.7,-70.7083], Tmin=(1392.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-511.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sF1s)) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
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
    label = 'F[C](F)C1(F)C=C[C]1F(2958)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {2,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-351.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,346,659,817,1284,2750,3150,900,1100,190,488,555,1236,1407,448.019,449.042,450.152,450.845,454.576,454.944,456.469],'cm^-1')),
        HinderedRotor(inertia=(0.000821279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16446,0.0655791,-7.03308e-05,3.89029e-08,-8.7669e-12,-42156.8,24.7383], Tmin=(100,'K'), Tmax=(1059.89,'K')), NASAPolynomial(coeffs=[12.0003,0.0246849,-1.24558e-05,2.49999e-09,-1.80442e-13,-44453.8,-28.1688], Tmin=(1059.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-351.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(CsCsF1sF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC1(F)OC2(F)C=CC12F(13762)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {7,S} {9,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-726.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09413,0.0681056,-6.35876e-05,2.94945e-08,-5.62755e-12,-87318.6,19.4449], Tmin=(100,'K'), Tmax=(1220.38,'K')), NASAPolynomial(coeffs=[13.0708,0.0288499,-1.53373e-05,3.13626e-09,-2.27949e-13,-90241.8,-40.721], Tmin=(1220.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-726.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(CsCCFO) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'OC(F)(F)C1(F)C=C=C1F(13765)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u0 p0 c0 {6,S} {10,D} {11,S}
9  C u0 p0 c0 {4,S} {6,S} {10,D}
10 C u0 p0 c0 {8,D} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-533.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0743101,0.096658,-0.000165905,1.51117e-07,-5.39083e-11,-64063.6,28.0593], Tmin=(100,'K'), Tmax=(805.157,'K')), NASAPolynomial(coeffs=[9.39367,0.0339309,-1.8438e-05,3.67223e-09,-2.5833e-13,-65031.8,-11.5746], Tmin=(805.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-533.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCddCF) + group(Cdd-CdsCds) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[O]C(F)(F)C1(F)C2[CH]C21F(13895)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,S} {10,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
10 C u1 p0 c0 {7,S} {8,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-403.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658189,0.0810682,-0.000109828,8.12553e-08,-2.47524e-11,-48465.4,26.4822], Tmin=(100,'K'), Tmax=(793.06,'K')), NASAPolynomial(coeffs=[10.3484,0.0321939,-1.73879e-05,3.54904e-09,-2.57038e-13,-50002.5,-18.0209], Tmin=(793.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(Cs-CsCsCsH) + group(CsCCCF) + group(Cs-CsCsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + polycyclic(s2_3_3_ane) + radical(O2sj(Cs-CsF1sF1s)) + radical(bicyclo[1.1.0]butane-secondary) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1C2[CH]C1(F)C(F)(F)O2(13838)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u1 p0 c0 {4,S} {6,S} {7,S}
10 C u1 p0 c0 {6,S} {7,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-507.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.87116,0.0602792,-4.69054e-05,1.68797e-08,-2.37748e-12,-60910.6,25.9051], Tmin=(100,'K'), Tmax=(1683.14,'K')), NASAPolynomial(coeffs=[19.0284,0.0171281,-8.44941e-06,1.64786e-09,-1.15055e-13,-67022.9,-71.147], Tmin=(1683.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-507.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + polycyclic(s3_4_5_ane) + radical(CsCsCsF1s) + radical(CCJCO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]1[CH]C2OC(F)(F)C12F(13849)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u1 p0 c0 {4,S} {6,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-443.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35821,0.0625864,-5.45047e-05,2.35847e-08,-4.23828e-12,-53259.4,24.5391], Tmin=(100,'K'), Tmax=(1275.13,'K')), NASAPolynomial(coeffs=[11.9679,0.0293045,-1.53539e-05,3.1159e-09,-2.25222e-13,-55965.1,-29.2255], Tmin=(1275.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(CsCFFO) + polycyclic(s2_4_4_ane) + radical(CsCsCsF1s) + radical(CCJCO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[CH]=CC(F)=C(F)C([O])(F)F(13896)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {11,S}
10 C u1 p0 c0 {9,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-429.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,323,467,575,827,1418,280,518,736,852,873,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.958056,'amu*angstrom^2'), symmetry=1, barrier=(22.0276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958686,'amu*angstrom^2'), symmetry=1, barrier=(22.0421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0863969,0.0883538,-0.000111046,6.86482e-08,-1.66659e-11,-51490.2,29.469], Tmin=(100,'K'), Tmax=(1007.4,'K')), NASAPolynomial(coeffs=[17.0941,0.0208224,-1.0492e-05,2.10425e-09,-1.51969e-13,-54916.9,-52.7084], Tmin=(1007.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-429.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(O2sj(Cs-F1sF1sCd)) + radical(Cds_P)"""),
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
    label = '[O]C(F)(F)C1=C(F)C=C1(13897)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {6,D} {9,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u0 p0 c0 {7,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-194.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,280,518,736,852,873,2750,3150,900,1100,180,180,180,180,180,704.521,1064.29,1065.93,1067.98,4000],'cm^-1')),
        HinderedRotor(inertia=(1.33076,'amu*angstrom^2'), symmetry=1, barrier=(30.5967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476098,0.0691765,-7.23823e-05,3.62334e-08,-6.99128e-12,-23310.7,26.6508], Tmin=(100,'K'), Tmax=(1275.05,'K')), NASAPolynomial(coeffs=[18.7874,0.0117319,-4.80329e-06,8.99516e-10,-6.336e-14,-27980.3,-66.1399], Tmin=(1275.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F)) + radical(O2sj(Cs-F1sF1sCd))"""),
)

species(
    label = '[O][C](F)F(2580)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-255.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53788e-05,4.83907e-09,3.86802e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58999,0.00108969,-4.53602e-07,1.24872e-10,-1.13713e-14,-32130.4,-16.3888], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
)

species(
    label = 'O=C(F)C1(F)C=C[C]1F(13898)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u1 p0 c0 {2,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  C u0 p0 c0 {3,S} {4,D} {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-498.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,346,659,817,1284,2750,3150,900,1100,486,617,768,1157,1926,524.501,524.52,524.527,524.532,524.561,524.614,524.708],'cm^-1')),
        HinderedRotor(inertia=(0.000613554,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46085,0.0555731,-4.71437e-05,1.96973e-08,-3.34428e-12,-59833.1,24.3717], Tmin=(100,'K'), Tmax=(1372.61,'K')), NASAPolynomial(coeffs=[12.7778,0.0225937,-1.11033e-05,2.19271e-09,-1.5606e-13,-62939.8,-33.81], Tmin=(1372.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-498.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(COCsFO) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = '[O]C(F)(F)C1(F)C=C=C1F(13899)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u0 p0 c0 {4,S} {6,S} {10,D}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {8,D} {9,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-273.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,145,326,398,834,1303,2950,1000,180,180,180,1776.21,1776.24],'cm^-1')),
        HinderedRotor(inertia=(0.265388,'amu*angstrom^2'), symmetry=1, barrier=(6.10179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.287292,0.0926888,-0.000164758,1.53471e-07,-5.53999e-11,-32805,27.8835], Tmin=(100,'K'), Tmax=(815.272,'K')), NASAPolynomial(coeffs=[8.39809,0.0329609,-1.81901e-05,3.62949e-09,-2.54996e-13,-33465.1,-5.52711], Tmin=(815.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCddCF) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(O2sj(Cs-CsF1sF1s)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
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
    label = 'O=C(F)[C]1C=C[C]1F(13900)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {8,D}
4  C u1 p0 c0 {5,S} {6,S} {8,S}
5  C u1 p0 c0 {1,S} {4,S} {7,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {5,S} {6,D} {9,S}
8  C u0 p0 c0 {2,S} {3,D} {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-173.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([346,659,817,1284,2750,3150,900,1100,611,648,830,1210,1753,700.257,700.261,700.262,700.274,700.279,700.282,700.285,700.286,700.287,700.29],'cm^-1')),
        HinderedRotor(inertia=(0.000343782,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0059,0.0318711,1.25904e-05,-3.88967e-08,1.65156e-11,-20812.5,21.557], Tmin=(100,'K'), Tmax=(1016.69,'K')), NASAPolynomial(coeffs=[12.8413,0.0173736,-7.52619e-06,1.51068e-09,-1.12882e-13,-24469.7,-38.0471], Tmin=(1016.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(CsCCFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(COCsFO) + ring(Cs(F)-Cs-Cd-Cd) + radical(C=CCJ(C)C=O) + radical(CsCdCsF1s)"""),
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
    label = '[O]C(F)(F)C1=C(F)C=[C]1(13463)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {9,S}
7  C u0 p0 c0 {3,S} {6,D} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (24.5023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,280,518,736,852,873,2950,1000,180,180,180,180,180,180,498.952,1166.74,2892.41],'cm^-1')),
        HinderedRotor(inertia=(0.0208929,'amu*angstrom^2'), symmetry=1, barrier=(20.1205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767926,0.0699759,-8.40722e-05,4.88786e-08,-1.10491e-11,3064.39,26.0763], Tmin=(100,'K'), Tmax=(1085.75,'K')), NASAPolynomial(coeffs=[15.9939,0.0138827,-6.57851e-06,1.29697e-09,-9.32604e-14,-241.982,-48.6332], Tmin=(1085.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.5023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(Cds-Cds(Cds-Cds)Cs) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F)) + radical(O2sj(Cs-F1sF1sCd)) + radical(cyclobutadiene-C1)"""),
)

species(
    label = '[O]C(F)(F)C1(F)C=[C]C1F(13505)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-388.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,164,312,561,654,898,1207,1299,3167,351,323,533,609,664,892,1120,1201,2950,1000,180,180,180,930.696,930.746],'cm^-1')),
        HinderedRotor(inertia=(0.00879764,'amu*angstrom^2'), symmetry=1, barrier=(5.4086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747538,0.0760113,-8.77886e-05,5.24208e-08,-1.26918e-11,-46567.2,26.9009], Tmin=(100,'K'), Tmax=(993.424,'K')), NASAPolynomial(coeffs=[12.9118,0.0270319,-1.38325e-05,2.7899e-09,-2.0182e-13,-48984,-31.7042], Tmin=(993.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-388.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sF1s)) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)(F)C1(F)[C]=CC1F(13901)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u0 p0 c0 {7,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-389.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,351,323,533,609,664,892,1120,1201,2950,1000,180,180,180,180,1224.79,1225.01],'cm^-1')),
        HinderedRotor(inertia=(0.122255,'amu*angstrom^2'), symmetry=1, barrier=(2.81088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64373,0.0804063,-0.000105977,7.53792e-08,-2.19048e-11,-46680.9,27.2288], Tmin=(100,'K'), Tmax=(832.683,'K')), NASAPolynomial(coeffs=[11.0663,0.0303411,-1.57928e-05,3.17869e-09,-2.28638e-13,-48416.7,-21.1461], Tmin=(832.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sF1s)) + radical(cyclobutene-vinyl) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'OC(F)(F)C1(F)[C]=C[C]1F(13902)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {11,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-511.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482303,0.0822868,-0.00010273,6.58524e-08,-1.69476e-11,-61349.1,28.0888], Tmin=(100,'K'), Tmax=(942.216,'K')), NASAPolynomial(coeffs=[13.6846,0.0262391,-1.35019e-05,2.71898e-09,-1.9625e-13,-63836.9,-34.8186], Tmin=(942.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-511.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(cyclobutene-vinyl) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'OC(F)(F)C1(F)C=[C][C]1F(13903)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {6,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-510.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441611,0.0796511,-9.09353e-05,5.15322e-08,-1.15836e-11,-61229.3,28.2751], Tmin=(100,'K'), Tmax=(1078.88,'K')), NASAPolynomial(coeffs=[15.916,0.0222792,-1.11693e-05,2.24282e-09,-1.62225e-13,-64568.3,-47.5546], Tmin=(1078.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cs-CsF1sH)(Cd-CsH)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FOC(F)(F)[C]1C=C[C]1F(13904)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {12,S}
10 C u0 p0 c0 {8,S} {9,D} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-290.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,346,659,817,1284,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314287,0.069977,-5.97666e-05,2.20265e-08,-2.58688e-12,-34842.3,29.2301], Tmin=(100,'K'), Tmax=(1184.66,'K')), NASAPolynomial(coeffs=[19.7116,0.0167439,-7.88938e-06,1.56984e-09,-1.13725e-13,-40298.5,-71.2692], Tmin=(1184.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFH) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CCJ(C)CO) + radical(CsCdCsF1s)"""),
)

species(
    label = '[O]C(F)(F)[C]1C=CC1(F)F(13905)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {11,S}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-530.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483163,0.0669991,-5.90582e-05,2.47832e-08,-4.05763e-12,-63728.7,28.4109], Tmin=(100,'K'), Tmax=(1477.12,'K')), NASAPolynomial(coeffs=[19.3387,0.0159386,-7.20652e-06,1.38091e-09,-9.68148e-14,-69299,-69.9114], Tmin=(1477.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-530.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(CsCCFF) + group(CsCFFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(O2sj(Cs-CsF1sF1s)) + radical(CCJ(C)CO)"""),
)

species(
    label = 'FO[C](F)C1(F)C=C[C]1F(13906)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u1 p0 c0 {2,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u0 p0 c0 {7,S} {8,D} {12,S}
10 C u1 p0 c0 {3,S} {5,S} {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-212.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,346,659,817,1284,2750,3150,900,1100,395,473,707,1436,180,449.211,449.277,449.553,449.556,449.792,449.918,449.943],'cm^-1')),
        HinderedRotor(inertia=(0.0152905,'amu*angstrom^2'), symmetry=1, barrier=(2.19765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348138,'amu*angstrom^2'), symmetry=1, barrier=(49.8457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.412499,0.0859648,-0.000115301,8.18717e-08,-2.36024e-11,-25420,28.3548], Tmin=(100,'K'), Tmax=(840.893,'K')), NASAPolynomial(coeffs=[12.0964,0.0303877,-1.61644e-05,3.27703e-09,-2.36556e-13,-27385,-25.9888], Tmin=(840.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(CsCCFH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O][C](F)C1(F)C=CC1(F)F(13907)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u0 p0 c0 {7,S} {8,D} {12,S}
10 C u1 p0 c0 {4,S} {5,S} {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-487.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,274,345,380,539,705,1166,1213,2750,3150,900,1100,395,473,707,1436,212.922,212.957,212.991,213.015,213.031,1451.46,1451.56,1451.57],'cm^-1')),
        HinderedRotor(inertia=(0.00708257,'amu*angstrom^2'), symmetry=1, barrier=(10.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.485727,0.0845316,-0.000121447,9.27761e-08,-2.77517e-11,-58492.7,27.8904], Tmin=(100,'K'), Tmax=(667.372,'K')), NASAPolynomial(coeffs=[10.7731,0.0301023,-1.53606e-05,3.0339e-09,-2.14659e-13,-60026.8,-18.7861], Tmin=(667.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-487.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = '[O]C(F)(F)[C](F)C1C=C1F(13488)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-324.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,253,525,597,667,842,1178,1324,212,367,445,1450,323,467,575,827,1418,180,180,180,180,180,1367.95,1368.09,1368.24],'cm^-1')),
        HinderedRotor(inertia=(0.00483811,'amu*angstrom^2'), symmetry=1, barrier=(6.42684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00483866,'amu*angstrom^2'), symmetry=1, barrier=(6.42845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493655,0.08386,-0.000120049,9.38554e-08,-2.98289e-11,-38849.5,31.3739], Tmin=(100,'K'), Tmax=(765.977,'K')), NASAPolynomial(coeffs=[10.6834,0.0306465,-1.58389e-05,3.15298e-09,-2.24415e-13,-40410.5,-15.0688], Tmin=(765.977,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-324.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cd-Cd-Cs(C-F)) + radical(O2sj(Cs-CsF1sF1s)) + radical(CsCsCsF1s)"""),
)

species(
    label = 'FC1=CC2OC(F)(F)C12F(13760)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u0 p0 c0 {4,S} {6,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-725.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803584,0.0659476,-5.70461e-05,2.31811e-08,-3.72315e-12,-87158.8,22.5304], Tmin=(100,'K'), Tmax=(1475.89,'K')), NASAPolynomial(coeffs=[17.7605,0.0199911,-1.03394e-05,2.0837e-09,-1.49517e-13,-92164.1,-65.8775], Tmin=(1475.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-725.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFFO) + group(CdCsCdF) + group(Cds-CdsCsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'FC1(F)OC2(F)[CH][CH]C12F(13879)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
9  C u1 p0 c0 {6,S} {10,S} {12,S}
10 C u1 p0 c0 {7,S} {9,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-452.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14775,0.0612276,-4.82925e-05,1.77904e-08,-2.61901e-12,-54347.9,25.4392], Tmin=(100,'K'), Tmax=(1571.8,'K')), NASAPolynomial(coeffs=[16.0964,0.0231856,-1.19883e-05,2.39227e-09,-1.6989e-13,-59047.2,-53.4394], Tmin=(1571.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCFFO) + polycyclic(s2_4_4_ane) + radical(cyclobutane) + radical(CCJCO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C(F)(F)C(F)=CC=[C]F(13908)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {9,S} {11,S}
9  C u0 p0 c0 {8,S} {10,D} {12,S}
10 C u1 p0 c0 {4,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-442.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,323,467,575,827,1418,2995,3025,975,1000,1300,1375,400,500,1630,1680,167,640,1190,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.866664,'amu*angstrom^2'), symmetry=1, barrier=(19.9263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869824,'amu*angstrom^2'), symmetry=1, barrier=(19.999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0374098,0.0907483,-0.000120472,7.89225e-08,-2.01608e-11,-53049.8,29.805], Tmin=(100,'K'), Tmax=(963.048,'K')), NASAPolynomial(coeffs=[17.199,0.0191565,-8.96297e-06,1.73e-09,-1.22062e-13,-56369.7,-52.7016], Tmin=(963.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-442.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(CdCFH) + radical(O2sj(Cs-F1sF1sCd)) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = '[O]C(F)(F)C1=C=C[CH]1(13909)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {8,D}
6  C u1 p0 c0 {5,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u0 p0 c0 {5,D} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (197.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.065,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60621,0.0429714,-1.6724e-05,-1.16634e-08,7.84959e-12,23860.4,23.4272], Tmin=(100,'K'), Tmax=(1024.38,'K')), NASAPolynomial(coeffs=[13.9547,0.0153443,-6.42167e-06,1.25487e-09,-9.21257e-14,20250.1,-41.7179], Tmin=(1024.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(O2sj(Cs-F1sF1sCd)) + radical(C=CCJC=C)"""),
)

species(
    label = '[O]C(F)(F)C1(F)[C]=C=C1(13497)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
7  C u0 p0 c0 {5,S} {9,D} {10,S}
8  C u1 p0 c0 {5,S} {9,D}
9  C u0 p0 c0 {7,D} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (139.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,2950,1000,180,180,180,180,180,1803.1,1803.24],'cm^-1')),
        HinderedRotor(inertia=(0.174205,'amu*angstrom^2'), symmetry=1, barrier=(4.00532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.618221,0.0873677,-0.000166376,1.61223e-07,-5.91183e-11,16934,25.438], Tmin=(100,'K'), Tmax=(837.171,'K')), NASAPolynomial(coeffs=[6.36676,0.0322273,-1.79933e-05,3.575e-09,-2.49281e-13,16941.3,4.51822], Tmin=(837.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(O2sj(Cs-CsF1sF1s)) + radical(Cds_S)"""),
)

species(
    label = '[O]C(F)(F)C1(F)C[C]=C1F(13910)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u0 p0 c0 {4,S} {6,S} {10,D}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-393.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,351,323,533,609,664,892,1120,1201,246,474,533,1155,180,180,180,986.8,999.349,1003.27,1462.36],'cm^-1')),
        HinderedRotor(inertia=(0.101263,'amu*angstrom^2'), symmetry=1, barrier=(2.32823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.834339,0.0730906,-7.95444e-05,4.41344e-08,-9.92094e-12,-47218.3,26.3197], Tmin=(100,'K'), Tmax=(1065.22,'K')), NASAPolynomial(coeffs=[13.3572,0.0260667,-1.33282e-05,2.69366e-09,-1.95204e-13,-49886.3,-34.8872], Tmin=(1065.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-393.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd(F)-Cd-Cs-Cs) + radical(O2sj(Cs-CsF1sF1s)) + radical(Cdj(Cs-CsHH)(Cd-CsF1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = '[O]C(F)(F)C1=C(F)[CH]C1F(13465)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {12,S}
10 C u0 p0 c0 {4,S} {8,D} {9,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-516.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,526,555,698,907,1200,1145,1227,2950,1000,271,519,563,612,1379,180,180,180,180,209.531,1162.83,1163.42],'cm^-1')),
        HinderedRotor(inertia=(0.0226252,'amu*angstrom^2'), symmetry=1, barrier=(21.7183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.268879,0.0769657,-8.00967e-05,4.01742e-08,-7.86739e-12,-62032.9,28.6556], Tmin=(100,'K'), Tmax=(1243.89,'K')), NASAPolynomial(coeffs=[18.8287,0.0172831,-8.1264e-06,1.6018e-09,-1.15072e-13,-66650.2,-64.9353], Tmin=(1243.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-516.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + group(Cds-CdsCsCs) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-F1sF1sCd)) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring)"""),
)

species(
    label = '[O][C](F)C1(F)C(F)=CC1F(13911)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {12,S}
10 C u1 p0 c0 {4,S} {5,S} {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-451.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,284,328,853,1146,1135,1297,3239,323,467,575,827,1418,2950,1000,395,473,707,1436,180,180,180,1607.62,1608.93],'cm^-1')),
        HinderedRotor(inertia=(0.327522,'amu*angstrom^2'), symmetry=1, barrier=(7.53038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.417526,0.0855828,-0.000130239,1.09301e-07,-3.70701e-11,-54195.3,29.6489], Tmin=(100,'K'), Tmax=(758.321,'K')), NASAPolynomial(coeffs=[10.1663,0.0309456,-1.58059e-05,3.10982e-09,-2.18921e-13,-55581.5,-14.0773], Tmin=(758.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-451.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCCF) + group(CsCCFH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FOC(F)(F)C1(F)[C]=C[CH]1(13912)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {6,S} {9,S} {11,S}
9  C u0 p0 c0 {8,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-174.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,2750,3150,900,1100,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635257,0.0770887,-8.53195e-05,4.76916e-08,-1.07334e-11,-20919,25.4267], Tmin=(100,'K'), Tmax=(1067.95,'K')), NASAPolynomial(coeffs=[14.3459,0.0257357,-1.31912e-05,2.6656e-09,-1.93163e-13,-23847.5,-41.6205], Tmin=(1067.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCCF) + group(Cs-(Cds-Cds)CsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(cyclobutene-allyl) + radical(cyclobutene-vinyl)"""),
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
    E0 = (-220.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (94.9861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (82.7691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (167.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (174.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-212.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-195.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (24.4695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-65.4272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-158.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-8.44509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (180.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (68.2517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-101.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-143.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (220.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (117.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (133.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (183.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (77.2272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (26.3205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-75.2484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-191.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (127.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-52.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (153.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (14.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (118.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-212.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-158.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-21.4507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (471.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (350.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (89.4471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-45.1812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (33.2714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (180.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['CF2O(49)', 'FC1=CC=C1F(304)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(F)(F)[C](F)C1(F)C=C1(13893)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['[O]C(F)(F)C1=CC(F)[C]1F(13894)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(303.62,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['[O]C(F)(F)C1C=C(F)[C]1F(13758)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(388.839,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'F[C](F)C1(F)C=C[C]1F(2958)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['FC1(F)OC2(F)C=CC12F(13762)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['OC(F)(F)C1(F)C=C=C1F(13765)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['[O]C(F)(F)C1(F)C2[CH]C21F(13895)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.00057e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['F[C]1C2[CH]C1(F)C(F)(F)O2(13838)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.17321e+11,'s^-1'), n=0.14784, Ea=(155.424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_beta;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['F[C]1[CH]C2OC(F)(F)C12F(13849)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.83633e+10,'s^-1'), n=0.345735, Ea=(62.6686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(F)=C(F)C([O])(F)F(13896)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(F)(F)C1=C(F)C=C1(13897)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(20.0362,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C](F)F(2580)', 'FC1=CC=C1F(304)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.89178e-07,'m^3/(mol*s)'), n=3.53001, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4944253016374622, var=1.7828810760479818, Tref=1000.0, N=135, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=C(F)C1(F)C=C[C]1F(13898)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(40.6319,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CF2O(49)', 'FC1=C(F)[CH][CH]1(923)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.33888e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(101.991,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]C(F)(F)C1(F)C=C=C1F(13899)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1090.54,'m^3/(mol*s)'), n=1.09894, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-5R!H-R_Ext-5R!H-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-5R!H-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C](F)F(2580)', 'FC1=C(F)[CH][CH]1(923)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(78)', 'O=C(F)[C]1C=C[C]1F(13900)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(33.2393,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[O]C(F)(F)C1=C(F)C=[C]1(13463)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(157.779,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(F)(F)C1(F)C=[C]C1F(13505)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(F)(F)C1(F)[C]=CC1F(13901)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.1e+10,'s^-1'), n=0.78, Ea=(132.633,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['OC(F)(F)C1(F)[C]=C[C]1F(13902)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['OC(F)(F)C1(F)C=[C][C]1F(13903)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.658e+09,'s^-1'), n=0.699, Ea=(29.5516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_doubleC] for rate rule [R5HJ_3;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['FOC(F)(F)[C]1C=C[C]1F(13904)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(135.843,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['[O]C(F)(F)[C]1C=CC1(F)F(13905)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(167.973,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FO[C](F)C1(F)C=C[C]1F(13906)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(83.3858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C](F)C1(F)C=CC1(F)F(13907)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(219.112,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(F)(F)[C](F)C1C=C1F(13488)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['FC1=CC2OC(F)(F)C12F(13760)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['FC1(F)OC2(F)[CH][CH]C12F(13879)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.83633e+10,'s^-1'), n=0.345735, Ea=(62.6686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(F)(F)C(F)=CC=[C]F(13908)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingle]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['F2(78)', '[O]C(F)(F)C1=C=C[CH]1(13909)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['HF(38)', '[O]C(F)(F)C1(F)[C]=C=C1(13497)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(209.062,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(F)(F)C1(F)C[C]=C1F(13910)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    products = ['[O]C(F)(F)C1=C(F)[CH]C1F(13465)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(175.67,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O][C](F)C1(F)C(F)=CC1F(13911)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(202.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction37',
    reactants = ['FOC(F)(F)C1(F)[C]=C[CH]1(13912)'],
    products = ['[O]C(F)(F)C1(F)C=C[C]1F(13756)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(72.7413,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #3817',
    isomers = [
        '[O]C(F)(F)C1(F)C=C[C]1F(13756)',
    ],
    reactants = [
        ('CF2O(49)', 'FC1=CC=C1F(304)'),
        ('CF2O(49)', 'FC1=C(F)[CH][CH]1(923)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3817',
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

