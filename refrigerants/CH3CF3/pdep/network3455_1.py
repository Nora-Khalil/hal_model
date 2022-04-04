species(
    label = 'F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {12,S} {14,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-573.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,212,367,445,1450,562,600,623,1070,1265,1685,370,180,1026.36,1026.38,1026.46,1026.5,1026.6],'cm^-1')),
        HinderedRotor(inertia=(0.17847,'amu*angstrom^2'), symmetry=1, barrier=(4.10338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178438,'amu*angstrom^2'), symmetry=1, barrier=(4.10263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.248578,0.091269,-0.000109832,7.22836e-08,-1.98811e-11,-68821,34.1605], Tmin=(100,'K'), Tmax=(867.757,'K')), NASAPolynomial(coeffs=[11.3411,0.040137,-2.14459e-05,4.38e-09,-3.18273e-13,-70746.2,-17.7814], Tmin=(867.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
10 C u1 p0 c0 {7,S} {8,S} {14,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-572.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([289,311,382,485,703,1397,215,315,519,588,595,1205,1248,164,312,561,654,898,1207,1299,3167,2950,1000,562,600,623,1070,1265,1685,370,180,1180.3,1185.88,1187.27],'cm^-1')),
        HinderedRotor(inertia=(0.175194,'amu*angstrom^2'), symmetry=1, barrier=(4.02805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174127,'amu*angstrom^2'), symmetry=1, barrier=(4.00352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.239617,0.102848,-0.000149338,1.21667e-07,-4.09483e-11,-68670,35.6444], Tmin=(100,'K'), Tmax=(720.22,'K')), NASAPolynomial(coeffs=[10.689,0.0421536,-2.29334e-05,4.66471e-09,-3.36115e-13,-70244.3,-13.4936], Tmin=(720.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCCF) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-F1sCsCs)(Cs-CsF1sF1s)(H)_ring) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,S} {12,D}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 C u1 p0 c0 {4,S} {9,S} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {9,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-723.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,350,440,435,1725,212,367,445,1450,234,589,736,816,1240,3237,182,240,577,636,1210,1413,180,1013.68,1013.74,1013.77,1013.77],'cm^-1')),
        HinderedRotor(inertia=(0.00310939,'amu*angstrom^2'), symmetry=1, barrier=(2.26737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985899,'amu*angstrom^2'), symmetry=1, barrier=(2.26678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00352308,0.0934212,-0.000109755,6.69291e-08,-1.65379e-11,-86869.3,36.0761], Tmin=(100,'K'), Tmax=(975.1,'K')), NASAPolynomial(coeffs=[14.7341,0.0329936,-1.67986e-05,3.37488e-09,-2.43504e-13,-89742.1,-34.6191], Tmin=(975.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-723.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(CdCFF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Csj(Cd-CsCd)(F1s)(H)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[C]=C(F)F(1205)',
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
    label = 'F[CH]C1[C](F)C1(F)F(10033)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-319.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,334,575,1197,1424,3202,282.23,715.634,715.834,716.386,716.536],'cm^-1')),
        HinderedRotor(inertia=(0.000334951,'amu*angstrom^2'), symmetry=1, barrier=(0.122425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39514,0.0592183,-6.43375e-05,3.56262e-08,-7.92326e-12,-38342.5,23.1859], Tmin=(100,'K'), Tmax=(1082.11,'K')), NASAPolynomial(coeffs=[12.0667,0.0197709,-9.65633e-06,1.93816e-09,-1.40322e-13,-40652.1,-29.1405], Tmin=(1082.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + radical(Csj(Cs-CsCsH)(F1s)(H)_1903_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FC(F)=C1C(F)C2C(F)(F)C12F(11826)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,S} {14,S}
11 C u0 p0 c0 {8,S} {10,S} {12,D}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-996.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25208,0.0919242,-0.000110391,7.19095e-08,-1.96063e-11,-119709,23.099], Tmin=(100,'K'), Tmax=(873.484,'K')), NASAPolynomial(coeffs=[11.4708,0.0405496,-2.21675e-05,4.57488e-09,-3.34403e-13,-121669,-29.5073], Tmin=(873.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-996.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCCCF) + group(CsCsCsFF) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC(F)=CC(F)C1=C(F)C1(F)F(11827)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 C u0 p0 c0 {8,S} {12,D} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-861.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.560176,0.110597,-0.000178086,1.57671e-07,-5.60573e-11,-103458,33.3111], Tmin=(100,'K'), Tmax=(766.887,'K')), NASAPolynomial(coeffs=[10.6996,0.0414093,-2.23036e-05,4.46547e-09,-3.16732e-13,-104878,-16.0175], Tmin=(766.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-861.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F)) + ring(Cs-Cd(C)-Cd(F))"""),
)

species(
    label = 'FC(F)=C=C(F)C1C(F)C1(F)F(11828)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {7,S} {8,S} {14,S}
10 C u0 p0 c0 {4,S} {7,S} {12,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {10,D} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-871.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.097867,0.0914289,-0.000100833,5.57435e-08,-1.22823e-11,-104653,32.1743], Tmin=(100,'K'), Tmax=(1097.11,'K')), NASAPolynomial(coeffs=[17.4155,0.0275765,-1.35328e-05,2.6949e-09,-1.94056e-13,-108495,-53.9405], Tmin=(1097.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-871.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'F[C](F)C(F)=CC(F)[C]=C(F)F(3745)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {9,D} {14,S}
9  C u0 p0 c0 {2,S} {8,D} {10,S}
10 C u1 p0 c0 {3,S} {4,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {7,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-716.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,271,519,563,612,1379,161,297,490,584,780,1358,562,600,623,1070,1265,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.439899,'amu*angstrom^2'), symmetry=1, barrier=(10.1142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438562,'amu*angstrom^2'), symmetry=1, barrier=(10.0834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9423,'amu*angstrom^2'), symmetry=1, barrier=(44.6573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.992195,0.12386,-0.000218481,2.00294e-07,-7.10425e-11,-86044.8,37.803], Tmin=(100,'K'), Tmax=(824.051,'K')), NASAPolynomial(coeffs=[10.7508,0.0415287,-2.25084e-05,4.44879e-09,-3.10542e-13,-87120.2,-11.3593], Tmin=(824.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + radical(Csj(Cd-F1sCd)(F1s)(F1s)) + radical(Cds_S)"""),
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
    label = 'FC(F)=[C]C(F)C1=C(F)C1(F)F(11829)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {8,S} {11,D}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-623.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,323,467,575,827,1418,562,600,623,1070,1265,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55428,'amu*angstrom^2'), symmetry=1, barrier=(35.7361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5549,'amu*angstrom^2'), symmetry=1, barrier=(35.7502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.658554,0.115996,-0.000205738,1.90577e-07,-6.84267e-11,-74851.9,34.1355], Tmin=(100,'K'), Tmax=(816.393,'K')), NASAPolynomial(coeffs=[9.85492,0.0402982,-2.22182e-05,4.42681e-09,-3.10668e-13,-75762.6,-9.51696], Tmin=(816.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-623.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCCFH) + group(Cds-CdsCsCs) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cd(C)-Cd(F)) + radical(Cds_S) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cds(F))"""),
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
    label = 'FC(F)=[C]C(F)C1C(F)=C1F(11830)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {3,S} {6,S} {8,D}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-342.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,164,312,561,654,898,1207,1299,3167,260,386,409,525,515,635,761,893,1354,1482,562,600,623,1070,1265,1685,370,180,180,1886.33,1888.25],'cm^-1')),
        HinderedRotor(inertia=(0.220244,'amu*angstrom^2'), symmetry=1, barrier=(5.06385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220182,'amu*angstrom^2'), symmetry=1, barrier=(5.06242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990277,0.0816409,-6.55844e-05,-6.94261e-08,1.1512e-10,-41063.8,27.3246], Tmin=(100,'K'), Tmax=(459.365,'K')), NASAPolynomial(coeffs=[9.27438,0.0387268,-2.08719e-05,4.16225e-09,-2.93662e-13,-42133.2,-9.5529], Tmin=(459.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-342.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'F[C]1C(C=C=C(F)F)C1(F)F(11831)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u0 p0 c0 {6,S} {11,D} {13,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {9,D} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-489.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,770.919,781.206,802.956,1174.3,2126.56],'cm^-1')),
        HinderedRotor(inertia=(0.207191,'amu*angstrom^2'), symmetry=1, barrier=(5.13339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.501627,0.0804768,-8.91937e-05,5.10929e-08,-1.18537e-11,-58800.4,32.3757], Tmin=(100,'K'), Tmax=(1035.53,'K')), NASAPolynomial(coeffs=[13.8981,0.0287303,-1.42385e-05,2.83805e-09,-2.04153e-13,-61574.9,-32.7226], Tmin=(1035.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-489.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCddFF) + group(Cdd-CdsCds) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'F[C]1C(C(F)=C=C(F)F)C1(F)F(11832)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {12,D}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u0 p0 c0 {10,D} {11,D}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-680.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,212,367,445,1450,145,326,398,834,1303,94,120,354,641,825,1294,540,610,2055,616.461,619.798,827.003,2367.2,2368.83],'cm^-1')),
        HinderedRotor(inertia=(0.00728016,'amu*angstrom^2'), symmetry=1, barrier=(28.965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (187.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00391676,0.0946087,-0.000125304,8.654e-08,-2.40556e-11,-81752.3,34.526], Tmin=(100,'K'), Tmax=(874.273,'K')), NASAPolynomial(coeffs=[13.9237,0.0309238,-1.60412e-05,3.22479e-09,-2.31922e-13,-84186.3,-30.759], Tmin=(874.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-680.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CdCddCF) + group(CdCddFF) + group(Cdd-CdsCds) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'FC#CC(F)C1[C](F)C1(F)F(11833)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  C u0 p0 c0 {7,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {10,S} {13,S}
9  C u1 p0 c0 {4,S} {6,S} {7,S}
10 C u0 p0 c0 {8,S} {11,T}
11 C u0 p0 c0 {5,S} {10,T}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-377.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,233,378,609,1068,1270,1314,3037,212,367,445,1450,2175,525,239,401,1367,180,771.066,771.091,771.094,771.109,771.111],'cm^-1')),
        HinderedRotor(inertia=(0.140947,'amu*angstrom^2'), symmetry=1, barrier=(3.24064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00768048,'amu*angstrom^2'), symmetry=1, barrier=(3.24073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164777,0.0912268,-0.000124721,9.27796e-08,-2.80635e-11,-45251.7,30.0627], Tmin=(100,'K'), Tmax=(803.53,'K')), NASAPolynomial(coeffs=[11.7622,0.0334916,-1.69377e-05,3.35001e-09,-2.38042e-13,-47115.4,-23.3506], Tmin=(803.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFH) + group(Ct-CtCs) + group(CtCF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'FC(F)=[C]C(F)[C]1C(F)=C1F(11834)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
7  C u1 p0 c0 {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u1 p0 c0 {6,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-210.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,206,336,431,607,515,611,528,696,1312,1446,562,600,623,1070,1265,1685,370,180,1409.43,1411.86],'cm^-1')),
        HinderedRotor(inertia=(0.18083,'amu*angstrom^2'), symmetry=1, barrier=(4.15764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184276,'amu*angstrom^2'), symmetry=1, barrier=(4.23686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.203476,0.0938247,-0.000157306,1.44905e-07,-5.26836e-11,-25155.7,28.139], Tmin=(100,'K'), Tmax=(798.494,'K')), NASAPolynomial(coeffs=[8.01368,0.0378043,-2.033e-05,4.04403e-09,-2.84857e-13,-25864.3,-4.41037], Tmin=(798.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(CsCCFH) + group(CdCsCdF) + group(CdCsCdF) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cd(F)-Cd) + radical(Allyl_T) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F)) + longDistanceInteraction_cyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'F[C]1[C](C=C=C(F)F)C1(F)F(11835)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u1 p0 c0 {6,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {11,D}
11 C u0 p0 c0 {9,D} {10,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-357.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,212,367,445,1450,3010,987.5,1337.5,450,1655,94,120,354,641,825,1294,540,610,2055,180,1056.39,1056.76,1058.19],'cm^-1')),
        HinderedRotor(inertia=(0.00219144,'amu*angstrom^2'), symmetry=1, barrier=(1.73995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (168.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726555,0.0756863,-8.5778e-05,5.08984e-08,-1.22639e-11,-42935,30.5158], Tmin=(100,'K'), Tmax=(998.416,'K')), NASAPolynomial(coeffs=[12.7225,0.0276272,-1.35765e-05,2.68859e-09,-1.92569e-13,-45330.4,-27.3391], Tmin=(998.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(Cds-CdsCsH) + group(CdCddFF) + group(Cdd-CdsCds) + ring(Cs-Cs(F)(F)-Cs) + radical(Allyl_T) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'FC(F)=[C]C(F)[C]1C(F)C1(F)F(3753)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {10,S} {12,S} {14,S}
10 C u1 p0 c0 {7,S} {8,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {12,D}
12 C u1 p0 c0 {9,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-612.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.118758,0.102037,-0.000166143,1.55184e-07,-5.79286e-11,-73528,36.6613], Tmin=(100,'K'), Tmax=(782.672,'K')), NASAPolynomial(coeffs=[7.05673,0.0469303,-2.52013e-05,5.03813e-09,-3.56866e-13,-74086.6,7.40906], Tmin=(782.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-612.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Cdj(Cs-F1sCsH)(Cd-F1sF1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C](F)C=C(F)C1[C](F)C1(F)F(11836)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {10,D} {12,S} {14,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-741.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.079193,0.0968757,-0.000122355,8.13455e-08,-2.19568e-11,-89051.2,35.7726], Tmin=(100,'K'), Tmax=(895.693,'K')), NASAPolynomial(coeffs=[13.8841,0.0345171,-1.79224e-05,3.61503e-09,-2.60835e-13,-91552.5,-30.0541], Tmin=(895.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-741.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Csj(Cd-CdH)(F1s)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)[C]=C(F)C1C(F)C1(F)F(11837)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {7,S} {8,S} {14,S}
10 C u0 p0 c0 {4,S} {7,S} {12,D}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {10,D} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-672.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0486045,0.0925133,-0.000102694,5.73177e-08,-1.28499e-11,-80790,34.2713], Tmin=(100,'K'), Tmax=(1073.21,'K')), NASAPolynomial(coeffs=[16.6969,0.0301005,-1.54608e-05,3.12968e-09,-2.27056e-13,-84384.3,-47.6994], Tmin=(1073.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-672.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cs) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsF1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1[C](C(F)C=C(F)F)C1(F)F(11838)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {13,S}
9  C u1 p0 c0 {7,S} {8,S} {10,S}
10 C u1 p0 c0 {4,S} {7,S} {9,S}
11 C u0 p0 c0 {8,S} {12,D} {14,S}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-650.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.470671,0.110191,-0.00018393,1.68045e-07,-6.05091e-11,-78138.7,36.0722], Tmin=(100,'K'), Tmax=(803.904,'K')), NASAPolynomial(coeffs=[9.0518,0.043385,-2.30318e-05,4.55716e-09,-3.19942e-13,-79042,-3.88569], Tmin=(803.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]1[C](F)C1C(F)C(F)=C(F)F(11839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {7,S} {11,S} {14,S}
9  C u1 p0 c0 {3,S} {7,S} {10,S}
10 C u1 p0 c0 {4,S} {7,S} {9,S}
11 C u0 p0 c0 {2,S} {8,S} {12,D}
12 C u0 p0 c0 {5,S} {6,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-572.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,284,328,853,1146,1135,1297,3239,165,259,333,401,397,493,1395,1505,323,467,575,827,1418,182,240,577,636,1210,1413,282.08,282.41,282.577,282.619,1549,1549.27],'cm^-1')),
        HinderedRotor(inertia=(7.02158e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900785,'amu*angstrom^2'), symmetry=1, barrier=(50.9412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.469261,0.106927,-0.000153686,1.18347e-07,-3.68612e-11,-68757.7,34.5403], Tmin=(100,'K'), Tmax=(782.201,'K')), NASAPolynomial(coeffs=[13.2255,0.036895,-1.9387e-05,3.8848e-09,-2.77584e-13,-70900.1,-28.1649], Tmin=(782.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(F1s)_ring) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(F1s)_ring) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)[C]=CC1C(F)(F)C1(F)F(11840)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {12,D} {14,S}
11 C u1 p0 c0 {5,S} {6,S} {12,S}
12 C u1 p0 c0 {10,D} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-677.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.112811,0.0934837,-0.000104239,5.82071e-08,-1.30148e-11,-81393,34.5086], Tmin=(100,'K'), Tmax=(1077.94,'K')), NASAPolynomial(coeffs=[17.161,0.0293844,-1.50423e-05,3.04198e-09,-2.20634e-13,-85117,-50.1237], Tmin=(1077.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-677.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs-Cs(F)(F)-Cs) + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Cdj(Cs-F1sF1sH)(Cd-CsH)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2)"""),
)

species(
    label = 'F[C](F)C(F)=CC1[C](F)C1(F)F(11841)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {14,S}
11 C u0 p0 c0 {4,S} {10,D} {12,S}
12 C u1 p0 c0 {5,S} {6,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-709.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.147328,0.0984287,-0.000126518,8.59167e-08,-2.36603e-11,-85232.9,36.9451], Tmin=(100,'K'), Tmax=(879.379,'K')), NASAPolynomial(coeffs=[13.8462,0.0347762,-1.79424e-05,3.60347e-09,-2.59197e-13,-87694,-28.7669], Tmin=(879.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-709.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(Cds-CdsCsH) + group(CdCsCdF) + ring(Cs-Cs(F)(F)-Cs) + radical(CsCsCsF1s) + radical(Csj(Cd-F1sCd)(F1s)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]=[C]C(F)C1C(F)(F)C1(F)F(11842)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {7,S} {11,S} {14,S}
11 C u1 p0 c0 {10,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
"""),
    E0 = (-558.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,186,258,302,356,410,480,483,561,545,633,1193,1235,1433,1517,164,312,561,654,898,1207,1299,3167,1685,370,167,640,1190,180,1071.41,1071.44,1071.48,1071.71],'cm^-1')),
        HinderedRotor(inertia=(0.216904,'amu*angstrom^2'), symmetry=1, barrier=(4.98704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216825,'amu*angstrom^2'), symmetry=1, barrier=(4.98524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.248248,0.104358,-0.000146967,1.05386e-07,-2.68585e-11,-67015.3,33.8651], Tmin=(100,'K'), Tmax=(610.41,'K')), NASAPolynomial(coeffs=[11.6952,0.0396188,-2.11199e-05,4.24409e-09,-3.0317e-13,-68725.4,-19.9239], Tmin=(610.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-558.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCCFH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Cdj(Cs-F1sCsH)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)2)"""),
)

species(
    label = 'F[C]=C(F)C(F)C1[C](F)C1(F)F(11843)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {12,S}
7  C u0 p0 c0 {8,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {11,S} {14,S}
10 C u1 p0 c0 {4,S} {7,S} {8,S}
11 C u0 p0 c0 {5,S} {9,S} {12,D}
12 C u1 p0 c0 {6,S} {11,D}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-547.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,284,328,853,1146,1135,1297,3239,212,367,445,1450,246,474,533,1155,167,640,1190,318.763,318.822,318.83,1187.54,2396.56,2396.59,2396.6],'cm^-1')),
        HinderedRotor(inertia=(0.54515,'amu*angstrom^2'), symmetry=1, barrier=(39.3183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545388,'amu*angstrom^2'), symmetry=1, barrier=(39.3184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (188.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.463285,0.106515,-0.000152007,1.15722e-07,-3.55629e-11,-65750.4,34.4389], Tmin=(100,'K'), Tmax=(793.019,'K')), NASAPolynomial(coeffs=[13.4967,0.0361025,-1.88248e-05,3.76203e-09,-2.68495e-13,-67964.5,-29.6728], Tmin=(793.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-547.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFF) + group(CsCsCsFH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sF1s)(F1s)_ring) + radical(Cdj(Cd-CsF1s)(F1s)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-190.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-31.7148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-95.6798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (259.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-181.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-90.4644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-126.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-138.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-103.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-13.8308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (146.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (38.5886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-78.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-61.3798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (119.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (78.9788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (79.0292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (40.7433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-16.0998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-7.57309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-49.1296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-48.1775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (51.2029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-11.2953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-62.5819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (14.8509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (4.22073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['FC1=CC1(F)F(974)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC(F)=[C]C(F)C1(F)[CH]C1(F)F(11625)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[CH]C(=C(F)F)C1[C](F)C1(F)F(11556)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=C(F)F(1205)', 'F[CH]C1[C](F)C1(F)F(10033)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['FC(F)=C1C(F)C2C(F)(F)C12F(11826)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['FC(F)=CC(F)C1=C(F)C1(F)F(11827)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.47101e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['FC(F)=C=C(F)C1C(F)C1(F)F(11828)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C](F)C(F)=CC(F)[C]=C(F)F(3745)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FC1=CC1(F)F(974)', 'F[CH][C]=C(F)F(2842)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(42.8173,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'FC(F)=[C]C(F)C1=C(F)C1(F)F(11829)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(45.45,'m^3/(mol*s)'), n=1.71, Ea=(14.8877,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_N-Sp-9R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R_Ext-5R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-6R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-7R!H-R_Ext-7R!H-R_5R!H-inRing_Ext-5R!H-R_N-Sp-9R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC(F)=[C]C(F)C1C(F)=C1F(11830)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(33.0164,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[C]1C(C=C=C(F)F)C1(F)F(11831)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(72.4976,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['FC=C=C(F)F(1325)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000502707,'m^3/(mol*s)'), n=2.87982, Ea=(6.29821,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'F[C]1C(C(F)=C=C(F)F)C1(F)F(11832)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(60.9,'m^3/(mol*s)'), n=1.64, Ea=(24.5742,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-1COS-R_6R!H-inRing',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-1COS-R_6R!H-inRing"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'FC#CC(F)C1[C](F)C1(F)F(11833)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(40.622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[CH][C]=C(F)F(2842)', 'F[C]1[CH]C1(F)F(11362)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'FC(F)=[C]C(F)[C]1C(F)=C1F(11834)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.28222,'m^3/(mol*s)'), n=1.29695, Ea=(187.235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'F[C]1[C](C=C=C(F)F)C1(F)F(11835)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(296.677,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['FC(F)=[C]C(F)[C]1C(F)C1(F)F(3753)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.93363e+09,'s^-1'), n=1.033, Ea=(174.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_Cs2] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[C](F)C=C(F)C1[C](F)C1(F)F(11836)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[C](F)[C]=C(F)C1C(F)C1(F)F(11837)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.95012e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[C]1[C](C(F)C=C(F)F)C1(F)F(11838)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C]1[C](F)C1C(F)C(F)=C(F)F(11839)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00110172,'s^-1'), n=4.50663, Ea=(241.043,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-4R!H-R',), comment="""Estimated from node R4F_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[C](F)[C]=CC1C(F)(F)C1(F)F(11840)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(178.859,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    products = ['F[C](F)C(F)=CC1[C](F)C1(F)F(11841)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(127.573,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['F[C]=[C]C(F)C1C(F)(F)C1(F)F(11842)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.59412e+07,'s^-1'), n=1.15307, Ea=(190.117,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[C]=C(F)C(F)C1[C](F)C1(F)F(11843)'],
    products = ['F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(169.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #3455',
    isomers = [
        'F[C]1C(C(F)[C]=C(F)F)C1(F)F(11626)',
    ],
    reactants = [
        ('FC1=CC1(F)F(974)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3455',
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

