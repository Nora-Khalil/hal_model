species(
    label = '[CH2]CC1[CH]C1F(9973)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {3,S} {11,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (257.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,904.312,904.328,904.331,904.36,904.373,904.376,904.377],'cm^-1')),
        HinderedRotor(inertia=(0.204656,'amu*angstrom^2'), symmetry=1, barrier=(4.70544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204655,'amu*angstrom^2'), symmetry=1, barrier=(4.70542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8556,0.0364024,8.9469e-06,-3.46266e-08,1.48588e-11,31007.3,22.9842], Tmin=(100,'K'), Tmax=(1012.63,'K')), NASAPolynomial(coeffs=[11.3243,0.0239122,-9.45376e-06,1.78217e-09,-1.27783e-13,27812.4,-29.1224], Tmin=(1012.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(RCCJ)"""),
)

species(
    label = 'C2H4(30)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'FC1C=C1(7087)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {2,S} {4,D} {6,S}
4 C u0 p0 c0 {2,S} {3,D} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (53.2336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,725.22,725.471,725.523,725.698],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2885.09,'J/mol'), sigma=(4.85707,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=450.64 K, Pc=57.13 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.11981,-0.0104554,0.000119533,-1.9991e-07,1.08069e-10,6403.99,8.0387], Tmin=(10,'K'), Tmax=(590.766,'K')), NASAPolynomial(coeffs=[1.59625,0.0247584,-1.59029e-05,4.86605e-09,-5.67736e-13,6385.84,16.208], Tmin=(590.766,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(53.2336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1CDC1""", comment="""Thermo library: CHOF_G4"""),
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
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.00015498,3.26298e-06,-2.40422e-09,5.69498e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.56,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C1[CH]C1F(9896)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
4  C u1 p0 c0 {2,S} {3,S} {8,S}
5  C u1 p0 c0 {2,S} {9,S} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (264.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,3000,3100,440,815,1455,1000,180,180,663.731,1462.56,1463.24,1463.99],'cm^-1')),
        HinderedRotor(inertia=(0.046134,'amu*angstrom^2'), symmetry=1, barrier=(70.1505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0808,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41677,0.0279919,1.78784e-06,-2.02853e-08,9.33473e-12,31822.7,18.31], Tmin=(100,'K'), Tmax=(994.049,'K')), NASAPolynomial(coeffs=[8.95228,0.0175158,-6.27965e-06,1.13771e-09,-8.01621e-14,29741.7,-17.1131], Tmin=(994.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Isobutyl)"""),
)

species(
    label = 'FC1C2CCC12(9975)',
    structure = adjacencyList("""1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-69.0538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47926,0.0165416,6.62138e-05,-9.2752e-08,3.48368e-11,-8235.47,13.2775], Tmin=(100,'K'), Tmax=(987.958,'K')), NASAPolynomial(coeffs=[11.2319,0.0236895,-9.29519e-06,1.83061e-09,-1.37442e-13,-12043.2,-39.3612], Tmin=(987.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.0538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'CCC1=CC1F(9976)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-9.72459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94994,0.0339044,1.54927e-05,-4.13268e-08,1.73093e-11,-1085.93,20.6046], Tmin=(100,'K'), Tmax=(998.359,'K')), NASAPolynomial(coeffs=[11.0749,0.0238312,-9.16761e-06,1.71411e-09,-1.22891e-13,-4227.9,-30.0138], Tmin=(998.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.72459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd(C)-Cd-Cs(F))"""),
)

species(
    label = 'C=CC1CC1F(9977)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-54.3316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22892,0.0253171,3.94836e-05,-6.50693e-08,2.52858e-11,-6458.73,19.2665], Tmin=(100,'K'), Tmax=(992.648,'K')), NASAPolynomial(coeffs=[10.8671,0.0240707,-9.34925e-06,1.78833e-09,-1.30853e-13,-9827.21,-30.6731], Tmin=(992.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.3316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs)"""),
)

species(
    label = 'CCC1C=C1F(9974)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {2,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (48.6152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59841,0.0468212,-2.9163e-05,8.97885e-09,-1.11781e-12,5938.44,18.0827], Tmin=(100,'K'), Tmax=(1831.77,'K')), NASAPolynomial(coeffs=[13.0703,0.0217708,-8.65019e-06,1.51342e-09,-9.89536e-14,1735.59,-44.2065], Tmin=(1831.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.6152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cd(F)-Cd)"""),
)

species(
    label = '[CH2]CC=C[CH]F(8670)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {6,S} {10,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  C u1 p0 c0 {1,S} {4,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (108.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,234,589,736,816,1240,3237,953.956],'cm^-1')),
        HinderedRotor(inertia=(0.183786,'amu*angstrom^2'), symmetry=1, barrier=(4.22559,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185638,'amu*angstrom^2'), symmetry=1, barrier=(4.26819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59444,'amu*angstrom^2'), symmetry=1, barrier=(36.6593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43311,0.0488627,-3.16365e-05,9.97279e-09,-1.25731e-12,13143.3,24.0305], Tmin=(100,'K'), Tmax=(1835.23,'K')), NASAPolynomial(coeffs=[14.6292,0.0201008,-8.12821e-06,1.43313e-09,-9.40081e-14,8299.75,-47.6453], Tmin=(1835.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(CsCFHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Csj(Cd-CdH)(F1s)(H))"""),
)

species(
    label = '[CH2][CH2](138)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,971.411,1436.01,2688.07,2689.99],'cm^-1')),
        HinderedRotor(inertia=(0.0117579,'amu*angstrom^2'), symmetry=1, barrier=(17.2429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40736,0.0100313,6.40915e-06,-1.41289e-08,5.92664e-12,38288.2,6.11704], Tmin=(100,'K'), Tmax=(954.263,'K')), NASAPolynomial(coeffs=[5.5225,0.00856171,-2.90742e-06,5.02351e-10,-3.4457e-14,37547.7,-5.75282], Tmin=(954.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,25474.2,-0.444973], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]CC1=CC1F(9978)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {2,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {10,S}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (195.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,1000,3000,3100,440,815,1455,1000,497.814,497.815,497.817,497.822,497.831,497.842],'cm^-1')),
        HinderedRotor(inertia=(0.000680224,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000680233,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98783,0.0346277,5.09029e-06,-2.91588e-08,1.29189e-11,23596.7,22.2553], Tmin=(100,'K'), Tmax=(1007.95,'K')), NASAPolynomial(coeffs=[10.8703,0.0215348,-8.39826e-06,1.57132e-09,-1.12269e-13,20680.6,-26.2508], Tmin=(1007.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(CsCCFH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cd(C)-Cd-Cs(F)) + radical(RCCJ)"""),
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
    label = '[CH2]CC1C=C1(9979)',
    structure = adjacencyList("""multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {1,S} {3,D} {10,S}
5  C u1 p0 c0 {2,S} {11,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (422.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.109,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97883,0.0380562,-1.36975e-05,-5.85434e-09,4.06334e-12,50898.9,16.9354], Tmin=(100,'K'), Tmax=(1074.66,'K')), NASAPolynomial(coeffs=[9.41598,0.0224765,-8.84344e-06,1.61344e-09,-1.11648e-13,48601.6,-22.7314], Tmin=(1074.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC1C=C1F(9980)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {10,S}
6  C u1 p0 c0 {3,S} {11,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (253.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,3000,3100,440,815,1455,1000,180,1305.99,1306.08,1306.29,1306.82,1724.32],'cm^-1')),
        HinderedRotor(inertia=(0.385322,'amu*angstrom^2'), symmetry=1, barrier=(8.8593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385225,'amu*angstrom^2'), symmetry=1, barrier=(8.85707,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91993,0.0445006,-3.01721e-05,1.0489e-08,-1.51952e-12,30608,18.6936], Tmin=(100,'K'), Tmax=(1548.14,'K')), NASAPolynomial(coeffs=[9.91879,0.0238335,-1.01476e-05,1.86593e-09,-1.27024e-13,28131.3,-23.3922], Tmin=(1548.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs-Cd(F)-Cd) + radical(RCCJ)"""),
)

species(
    label = 'FC1[CH][CH]1(9903)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {4,S} {6,S}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (308.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,2750,3150,900,1100,1800.38,1800.4,1800.51,1800.62],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0542,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60777,0.00508457,2.30875e-05,-2.54937e-08,7.79572e-12,37164.6,14.3883], Tmin=(100,'K'), Tmax=(1144.06,'K')), NASAPolynomial(coeffs=[3.87844,0.0162025,-7.30706e-06,1.43501e-09,-1.0285e-13,36313.1,9.59549], Tmin=(1144.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = 'C=CC1[CH]C1F(9918)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {2,S} {3,S} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (171.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,259,529,569,1128,1321,1390,3140,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,400.667,1417.76,1417.77,1417.81,1417.86,1417.87],'cm^-1')),
        HinderedRotor(inertia=(0.932364,'amu*angstrom^2'), symmetry=1, barrier=(21.4369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39065,0.0258992,2.14448e-05,-3.93223e-08,1.47315e-11,20702.5,21.0358], Tmin=(100,'K'), Tmax=(1046.32,'K')), NASAPolynomial(coeffs=[8.4987,0.0254294,-1.06833e-05,2.04781e-09,-1.46747e-13,18171.9,-14.6937], Tmin=(1046.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cs(F)-Cs-Cs) + radical(cyclopropane)"""),
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
    label = '[CH2]C[C]1C=C1(9981)',
    structure = adjacencyList("""multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u1 p0 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {4,D} {8,S}
4  C u0 p0 c0 {2,S} {3,D} {9,S}
5  C u1 p0 c0 {1,S} {10,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (554.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,3150,900,1100,3000,3100,440,815,1455,1000,180,180,1094.67,1094.68,1094.69,1094.89,1094.96,1094.98],'cm^-1')),
        HinderedRotor(inertia=(0.137081,'amu*angstrom^2'), symmetry=1, barrier=(3.15177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137079,'amu*angstrom^2'), symmetry=1, barrier=(3.15172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26723,0.0325492,-7.9236e-06,-8.90727e-09,4.79132e-12,66761.5,14.1526], Tmin=(100,'K'), Tmax=(1049.63,'K')), NASAPolynomial(coeffs=[8.05294,0.0216725,-8.34556e-06,1.50128e-09,-1.03065e-13,64931.5,-16.9719], Tmin=(1049.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Allyl_T) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]1CC1F(8701)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
4  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {3,S} {4,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (204.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23717,0.037349,-1.32939e-05,-2.12103e-10,6.5067e-13,24610.3,22.0509], Tmin=(100,'K'), Tmax=(1676.58,'K')), NASAPolynomial(coeffs=[10.292,0.026627,-1.13017e-05,2.018e-09,-1.32526e-13,20715.4,-24.532], Tmin=(1676.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Tertalkyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC1C[C]1F(9982)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (248.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03216,0.0314189,2.2513e-05,-4.82267e-08,1.96165e-11,29953.1,22.0307], Tmin=(100,'K'), Tmax=(996.17,'K')), NASAPolynomial(coeffs=[10.9572,0.0240573,-9.27961e-06,1.74479e-09,-1.25736e-13,26762.1,-28.0847], Tmin=(996.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C1[CH]C1F(9983)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
4  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {3,S} {13,S}
6  C u1 p0 c0 {2,S} {4,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (246.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00108,0.0345022,8.04632e-06,-2.90946e-08,1.17359e-11,29713.3,23.2344], Tmin=(100,'K'), Tmax=(1059.29,'K')), NASAPolynomial(coeffs=[9.80386,0.0265326,-1.11057e-05,2.11459e-09,-1.50605e-13,26854.3,-20.5512], Tmin=(1059.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring) + radical(Cs_S)"""),
)

species(
    label = '[CH2][CH]C1CC1F(9984)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
5  C u1 p0 c0 {2,S} {6,S} {11,S}
6  C u1 p0 c0 {5,S} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (216.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13001,0.0330154,6.70928e-06,-2.3941e-08,9.06122e-12,26125.6,23.4624], Tmin=(100,'K'), Tmax=(1116.34,'K')), NASAPolynomial(coeffs=[8.53906,0.0286291,-1.23604e-05,2.35522e-09,-1.66333e-13,23537.1,-13.3478], Tmin=(1116.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'CC[C]1[CH]C1F(9985)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {2,S} {3,S} {6,S}
6  C u1 p0 c0 {3,S} {5,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (237.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99152,0.0417036,-2.10314e-05,4.58855e-09,-3.48886e-13,28610.9,22.5632], Tmin=(100,'K'), Tmax=(1968.48,'K')), NASAPolynomial(coeffs=[14.2005,0.0211559,-8.6209e-06,1.48518e-09,-9.4417e-14,22978.7,-46.7038], Tmin=(1968.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(F)-Cs-Cs) + radical(Tertalkyl) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = 'CCC1[CH][C]1F(9986)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  C u1 p0 c0 {1,S} {2,S} {6,S}
6  C u1 p0 c0 {2,S} {5,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (278.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,212,367,445,1450,666.259,666.667,666.688,666.719,666.753,666.858,666.894],'cm^-1')),
        HinderedRotor(inertia=(0.00988101,'amu*angstrom^2'), symmetry=1, barrier=(3.11813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00988689,'amu*angstrom^2'), symmetry=1, barrier=(3.11825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90484,0.0328364,2.43432e-05,-5.44062e-08,2.29062e-11,33540.8,21.8001], Tmin=(100,'K'), Tmax=(975.11,'K')), NASAPolynomial(coeffs=[12.5368,0.0214618,-7.75169e-06,1.4421e-09,-1.05012e-13,29934.6,-37.0842], Tmin=(975.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(CsCsCsFH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cs(C)-Cs-Cs(F)) + radical(Csj(Cs-CsCsH)(Cs-CsHH)(F1s)_ring) + radical(Csj(Cs-CsCsH)(Cs-CsF1sH)(H)_ring)"""),
)

species(
    label = 'FCCC1[CH][CH]1(9987)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
5  C u1 p0 c0 {2,S} {6,S} {12,S}
6  C u1 p0 c0 {2,S} {5,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (248.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10097,0.0343256,1.41615e-06,-1.84044e-08,7.23477e-12,29957.4,23.579], Tmin=(100,'K'), Tmax=(1130.63,'K')), NASAPolynomial(coeffs=[8.47964,0.0280858,-1.19667e-05,2.25908e-09,-1.58462e-13,27471.5,-12.5925], Tmin=(1130.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(CsCsFHH) + ring(Cyclopropane) + radical(cyclopropane) + radical(cyclopropane)"""),
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
    E0 = (49.1169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (437.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (57.4012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (112.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (112.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (74.0902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (126.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (165.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (205.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (335.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (258.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (151.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (175.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (419.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (319.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (203.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (210.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (196.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (206.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (143.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (171.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (272.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['C2H4(30)', 'FC1C=C1(7087)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(18)', '[CH2]C1[CH]C1F(9896)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['FC1C2CCC12(9975)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['CCC1=CC1F(9976)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['C=CC1CC1F(9977)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['CCC1C=C1F(9974)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CC=C[CH]F(8670)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH2](138)', 'FC1C=C1(7087)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.20532e-12,'m^3/(mol*s)'), n=5.04048, Ea=(2.53045,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5477765617574784, var=1.0081445727216138, Tref=1000.0, N=71, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_5R!H-inRing
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[CH2]CC1=CC1F(9978)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(6.45698,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2]CC1C=C1(9979)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.15e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.8449,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[CH2]CC1C=C1F(9980)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(1.02232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C2H4(30)', 'FC1[CH][CH]1(9903)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00201083,'m^3/(mol*s)'), n=2.87982, Ea=(8.52254,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R!H-inRing_Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'C=CC1[CH]C1F(9918)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.42983,'m^3/(mol*s)'), n=1.92039, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.021512250841106063, var=1.4639131859069563, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH2](138)', 'FC1[CH][CH]1(9903)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2e+08,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', '[CH2]C[C]1C=C1(9981)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(562.232,'m^3/(mol*s)'), n=1.03051, Ea=(253.861,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['[CH2]C[C]1CC1F(8701)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2] for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['[CH2]CC1C[C]1F(9982)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.41e+08,'s^-1'), n=1.52, Ea=(161.502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy3;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['C[CH]C1[CH]C1F(9983)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['[CH2][CH]C1CC1F(9984)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.5088e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['CC[C]1[CH]C1F(9985)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4313e+06,'s^-1'), n=1.72764, Ea=(94.0044,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCC1[CH][C]1F(9986)'],
    products = ['[CH2]CC1[CH]C1F(9973)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.52029e+07,'s^-1'), n=1.31708, Ea=(101.589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CC1[CH]C1F(9973)'],
    products = ['FCCC1[CH][CH]1(9987)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.56499e-06,'s^-1'), n=5.16802, Ea=(223.544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.010618623676007603, var=1.2440368903510177, Tref=1000.0, N=3, data_mean=0.0, correlation='R4F',), comment="""Estimated from node R4F"""),
)

network(
    label = 'PDepNetwork #3455',
    isomers = [
        '[CH2]CC1[CH]C1F(9973)',
    ],
    reactants = [
        ('C2H4(30)', 'FC1C=C1(7087)'),
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

