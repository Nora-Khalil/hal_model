species(
    label = '[O]C(C(=O)F)C([O])(F)C=O(1906)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u1 p2 c0 {7,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-649.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,486,617,768,1157,1926,300.669,301.046,301.14,800.049],'cm^-1')),
        HinderedRotor(inertia=(0.270305,'amu*angstrom^2'), symmetry=1, barrier=(17.4214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111922,'amu*angstrom^2'), symmetry=1, barrier=(7.19451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271496,'amu*angstrom^2'), symmetry=1, barrier=(17.421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.414711,0.0850207,-0.000118863,8.68011e-08,-2.54071e-11,-77976.2,34.7369], Tmin=(100,'K'), Tmax=(833.439,'K')), NASAPolynomial(coeffs=[12.5773,0.0266464,-1.37999e-05,2.75971e-09,-1.97375e-13,-80003.5,-21.7244], Tmin=(833.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-649.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(COCsFO) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = 'O=CC(=O)F(335)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C([O])C(=O)F(1893)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {6,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
8  C u0 p0 c0 {2,S} {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-542.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,486,617,768,1157,1926,180,1161.51,3799.88],'cm^-1')),
        HinderedRotor(inertia=(0.384624,'amu*angstrom^2'), symmetry=1, barrier=(8.84326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82972,'amu*angstrom^2'), symmetry=1, barrier=(19.0769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4361.52,'J/mol'), sigma=(6.57568,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=681.26 K, Pc=34.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53589,0.0556679,-6.17899e-05,3.52064e-08,-8.05544e-12,-65151.4,28.5373], Tmin=(100,'K'), Tmax=(1054.95,'K')), NASAPolynomial(coeffs=[11.3587,0.0184229,-8.8318e-06,1.7397e-09,-1.24475e-13,-67223.9,-19.3775], Tmin=(1054.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-542.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(COCsFO) + radical(C=OCOJ) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]C(F)C([O])(F)C=O(1894)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {6,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8  C u0 p0 c0 {5,D} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-506.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.990669,'amu*angstrom^2'), symmetry=1, barrier=(22.7774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992998,'amu*angstrom^2'), symmetry=1, barrier=(22.831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4349.1,'J/mol'), sigma=(6.772,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=679.32 K, Pc=31.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.068,0.0694556,-0.000101351,7.84664e-08,-2.43424e-11,-60868.4,25.7348], Tmin=(100,'K'), Tmax=(788.367,'K')), NASAPolynomial(coeffs=[10.4096,0.0220562,-1.11615e-05,2.19598e-09,-1.55055e-13,-62341.3,-17.1113], Tmin=(788.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]C(=COF)C([O])C(=O)F(2425)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {9,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {12,S}
10 C u0 p0 c0 {1,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-408.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,430.062,430.069,430.101,430.219,4000],'cm^-1')),
        HinderedRotor(inertia=(0.122296,'amu*angstrom^2'), symmetry=1, barrier=(16.0515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122293,'amu*angstrom^2'), symmetry=1, barrier=(16.0503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122292,'amu*angstrom^2'), symmetry=1, barrier=(16.051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0547233,0.0881011,-0.000114724,7.21733e-08,-1.74963e-11,-48964.6,36.4022], Tmin=(100,'K'), Tmax=(1019.7,'K')), NASAPolynomial(coeffs=[18.7753,0.0142352,-6.06427e-06,1.13225e-09,-7.89493e-14,-52804.8,-54.8092], Tmin=(1019.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-408.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(COCsFO) + radical(C=OCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(O[CH]C(=O)F)C(=O)F(1907)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
8  C u1 p0 c0 {3,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-704.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,486,617,768,1157,1926,611,648,830,1210,1753,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4473.71,'J/mol'), sigma=(6.57407,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=698.78 K, Pc=35.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0610863,0.093451,-0.000143074,1.12695e-07,-3.52087e-11,-84656.2,35.4887], Tmin=(100,'K'), Tmax=(785.424,'K')), NASAPolynomial(coeffs=[13.3856,0.0255897,-1.34673e-05,2.68106e-09,-1.90105e-13,-86749.2,-25.5755], Tmin=(785.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(COCsFO) + radical(C=OCOJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C(F)(C=O)O[C](F)C=O(1911)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u1 p0 c0 {2,S} {3,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-689.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,280,501,1494,1531,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,240.055,240.059,240.059,240.063],'cm^-1')),
        HinderedRotor(inertia=(1.03433,'amu*angstrom^2'), symmetry=1, barrier=(42.299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227264,'amu*angstrom^2'), symmetry=1, barrier=(9.29373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03442,'amu*angstrom^2'), symmetry=1, barrier=(42.299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03447,'amu*angstrom^2'), symmetry=1, barrier=(42.299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4426.76,'J/mol'), sigma=(6.80439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.45 K, Pc=31.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0627585,0.097645,-0.000168627,1.54588e-07,-5.50738e-11,-82777.9,31.5212], Tmin=(100,'K'), Tmax=(823.11,'K')), NASAPolynomial(coeffs=[8.71625,0.0355371,-1.88967e-05,3.71443e-09,-2.58846e-13,-83523,-4.41538], Tmin=(823.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-689.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C=C(F)C([O])C(=O)F(2426)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {8,D}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {2,S} {4,D} {6,S}
9  C u0 p0 c0 {5,S} {7,D} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-529.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,323,467,575,827,1418,486,617,768,1157,1926,3010,987.5,1337.5,450,1655,403.649,403.678,403.763,403.861],'cm^-1')),
        HinderedRotor(inertia=(0.00103369,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100977,'amu*angstrom^2'), symmetry=1, barrier=(11.6751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9364,0.067403,-7.44811e-05,4.07229e-08,-8.81227e-12,-63589.7,30.9355], Tmin=(100,'K'), Tmax=(1120.54,'K')), NASAPolynomial(coeffs=[14.6234,0.0185442,-9.07659e-06,1.81034e-09,-1.30603e-13,-66657,-36.6543], Tmin=(1120.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-529.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COCsFO) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(F)([CH]C(=O)F)C=O(2427)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,D} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-523.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,611,648,830,1210,1753,279.068,282.622,287.613],'cm^-1')),
        HinderedRotor(inertia=(0.136373,'amu*angstrom^2'), symmetry=1, barrier=(7.56151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130773,'amu*angstrom^2'), symmetry=1, barrier=(7.5005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664308,'amu*angstrom^2'), symmetry=1, barrier=(36.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613411,0.0824927,-0.000137042,1.21754e-07,-4.25586e-11,-62825.3,30.9074], Tmin=(100,'K'), Tmax=(808.715,'K')), NASAPolynomial(coeffs=[9.10351,0.0285971,-1.49999e-05,2.94907e-09,-2.06129e-13,-63809.2,-5.84342], Tmin=(808.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-523.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(COCsFO) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'O=CC1(F)OOC1C(=O)F(1916)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-744.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689393,0.0629903,-5.27952e-05,2.07242e-08,-3.17509e-12,-89439.7,31.648], Tmin=(100,'K'), Tmax=(1567.13,'K')), NASAPolynomial(coeffs=[19.2164,0.015701,-7.53134e-06,1.46862e-09,-1.03286e-13,-95246.6,-66.0574], Tmin=(1567.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-744.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(COCsFO) + ring(Cs-Cs(F)-O2s-O2s)"""),
)

species(
    label = 'O=CC(O)(F)C(=O)C(=O)F(2428)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {12,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,D} {7,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-1047.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0239656,0.0959809,-0.000148902,1.19908e-07,-3.84627e-11,-125834,30.3311], Tmin=(100,'K'), Tmax=(764.349,'K')), NASAPolynomial(coeffs=[13.0245,0.0276859,-1.48573e-05,2.97709e-09,-2.11733e-13,-127828,-29.1116], Tmin=(764.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1047.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(COCFO)"""),
)

species(
    label = '[O]C(C(=O)F)C1(F)[CH]OO1(2429)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u1 p0 c0 {4,S} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-392.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.937166,0.0709705,-7.68107e-05,4.28621e-08,-9.73644e-12,-47080.4,32.2366], Tmin=(100,'K'), Tmax=(1052.41,'K')), NASAPolynomial(coeffs=[12.6419,0.026483,-1.34026e-05,2.69507e-09,-1.94741e-13,-49544,-24.8299], Tmin=(1052.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-392.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(COCsFO) + ring(Cs-Cs(F)-O2s-O2s) + radical(C=OCOJ) + radical(CCsJOO)"""),
)

species(
    label = '[O]C(F)(C=O)C1OO[C]1F(2430)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
9  C u1 p0 c0 {2,S} {4,S} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-346.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279898,0.0890225,-0.000138893,1.17571e-07,-3.99966e-11,-41519.4,32.085], Tmin=(100,'K'), Tmax=(761.128,'K')), NASAPolynomial(coeffs=[10.6917,0.030446,-1.58477e-05,3.13508e-09,-2.21105e-13,-42992.6,-14.5696], Tmin=(761.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + group(Cds-OdCsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(C=OCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)[CH]OOC1C(=O)F(2305)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-490.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.219195,0.0748613,-7.56095e-05,3.73765e-08,-6.93943e-12,-58883.2,32.2188], Tmin=(100,'K'), Tmax=(1513.88,'K')), NASAPolynomial(coeffs=[19.839,0.0110305,-1.63064e-06,7.18393e-11,1.55006e-15,-63715,-68.7684], Tmin=(1513.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-490.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(COCsFO) + ring(12dioxolane) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOOC)"""),
)

species(
    label = '[O]C1[C](F)OOC1(F)C=O(2254)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {5,S} {7,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {4,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-458.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.159016,0.0799685,-8.65436e-05,4.61266e-08,-9.38758e-12,-55041.2,28.8376], Tmin=(100,'K'), Tmax=(1299.61,'K')), NASAPolynomial(coeffs=[20.0347,0.0127976,-3.22427e-06,4.15126e-10,-2.28101e-14,-59866.2,-72.2473], Tmin=(1299.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-458.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCFHO) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1OC1(F)C([O])C(=O)F(2431)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {5,S} {7,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-587.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26696,0.0669295,-7.01224e-05,4e-08,-9.71923e-12,-70554.4,30.9967], Tmin=(100,'K'), Tmax=(962.636,'K')), NASAPolynomial(coeffs=[9.38943,0.0331784,-1.75307e-05,3.57791e-09,-2.60272e-13,-72118.2,-7.88018], Tmin=(962.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-587.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(COCsFO) + ring(Cs(O2)-O2s-Cs(F)) + radical(C=OCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(F)(C=O)C1OC1([O])F(2432)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-527.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619623,0.0805207,-0.000104324,7.22967e-08,-2.04012e-11,-63282.6,29.1991], Tmin=(100,'K'), Tmax=(857.703,'K')), NASAPolynomial(coeffs=[11.5423,0.0295821,-1.52408e-05,3.05606e-09,-2.19498e-13,-65156.3,-21.82], Tmin=(857.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFOO) + group(Cds-OdCsH) + ring(Cs(F)(O2)-O2s-Cs) + radical(C=OCOJ) + radical(O2sj(Cs-F1sO2sCs))"""),
)

species(
    label = '[O]C1OC(C(=O)F)C1([O])F(2408)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {5,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-606.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38917,0.0563575,-4.17084e-05,1.4158e-08,-1.91147e-12,-72825.3,31.8794], Tmin=(100,'K'), Tmax=(1697.19,'K')), NASAPolynomial(coeffs=[15.929,0.0220896,-1.14218e-05,2.26123e-09,-1.59049e-13,-77760.7,-45.9577], Tmin=(1697.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-606.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(COCsFO) + ring(Cs-O2s-Cs-Cs(F)) + radical(O2sj(Cs-F1sCsCs)) + radical(CCOJ)"""),
)

species(
    label = '[O]C1C([O])(F)OC1(F)C=O(2324)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-584.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0507,0.0693668,-7.02483e-05,3.65035e-08,-7.81773e-12,-70157.4,27.0603], Tmin=(100,'K'), Tmax=(1101.07,'K')), NASAPolynomial(coeffs=[12.2161,0.0288047,-1.49899e-05,3.04606e-09,-2.2113e-13,-72616.2,-27.8814], Tmin=(1101.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-584.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCFOO) + group(Cds-OdCsH) + ring(O2s-Cs-Cs-Cs(F)) + radical(CC(C)OJ) + radical(O2sj(Cs-F1sO2sCs))"""),
)

species(
    label = '[O]C(F)(C=O)C(=O)[C](O)F(2433)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {9,S} {12,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
8  C u0 p0 c0 {5,D} {7,S} {9,S}
9  C u1 p0 c0 {2,S} {3,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-721.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,280,501,1494,1531,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.366101,'amu*angstrom^2'), symmetry=1, barrier=(8.41739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18121,'amu*angstrom^2'), symmetry=1, barrier=(50.1504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02665,'amu*angstrom^2'), symmetry=1, barrier=(23.6047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1732,'amu*angstrom^2'), symmetry=1, barrier=(49.9661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139991,0.102073,-0.000178363,1.60456e-07,-5.56073e-11,-86663.8,31.9745], Tmin=(100,'K'), Tmax=(840.86,'K')), NASAPolynomial(coeffs=[10.2627,0.0320318,-1.67486e-05,3.24944e-09,-2.23945e-13,-87686.6,-12.0882], Tmin=(840.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-721.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCOF1sO2s)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(C(=O)F)C(=O)C=O(2434)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {3,D} {6,S} {9,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-573.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,486,617,768,1157,1926,2782.5,750,1395,475,1775,1000,286.431,1685.66,1686.39],'cm^-1')),
        HinderedRotor(inertia=(0.161061,'amu*angstrom^2'), symmetry=1, barrier=(9.53666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165047,'amu*angstrom^2'), symmetry=1, barrier=(9.52596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.674525,'amu*angstrom^2'), symmetry=1, barrier=(39.4591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05048,0.0707889,-0.000102337,8.41857e-08,-2.85255e-11,-68847.1,30.2694], Tmin=(100,'K'), Tmax=(716.902,'K')), NASAPolynomial(coeffs=[8.53317,0.0290441,-1.50039e-05,2.98308e-09,-2.1199e-13,-69920.1,-3.34083], Tmin=(716.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCsFO) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(509)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.101,381.695],'cm^-1')),
        HinderedRotor(inertia=(0.482775,'amu*angstrom^2'), symmetry=1, barrier=(49.7784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80088,0.0290439,-3.02832e-05,1.66615e-08,-3.81631e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.22,'K')), NASAPolynomial(coeffs=[6.95396,0.0128879,-6.71487e-06,1.38086e-09,-1.01081e-13,-26608.7,-6.32755], Tmin=(1028.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
)

species(
    label = 'HCO(15)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C(C(=O)F)C(=O)F(2235)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {2,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-719.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,430,542,558,676,687,849,1103,1211,1906,1946,180,180,1679.68],'cm^-1')),
        HinderedRotor(inertia=(0.0991169,'amu*angstrom^2'), symmetry=1, barrier=(2.27889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995248,'amu*angstrom^2'), symmetry=1, barrier=(2.28827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74736,0.0539537,-7.42851e-05,5.72839e-08,-1.81603e-11,-86430.5,25.5702], Tmin=(100,'K'), Tmax=(764.883,'K')), NASAPolynomial(coeffs=[7.90924,0.0217294,-1.10896e-05,2.20241e-09,-1.5679e-13,-87373.1,-2.50564], Tmin=(764.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(COCsFO) + group(COCsFO) + radical(C=OCOJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)(C=O)C=O(2195)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 C u0 p0 c0 {4,D} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-432.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(1.14372,'amu*angstrom^2'), symmetry=1, barrier=(26.2963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14945,'amu*angstrom^2'), symmetry=1, barrier=(26.4281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80007,0.0544787,-7.10848e-05,4.29746e-08,-5.76321e-12,-51904.3,20.4612], Tmin=(100,'K'), Tmax=(585.069,'K')), NASAPolynomial(coeffs=[8.00131,0.0215576,-1.0975e-05,2.16307e-09,-1.52774e-13,-52792.1,-7.51811], Tmin=(585.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-432.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
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
    label = '[O]C(F)(C=O)C(=O)C(=O)F(2435)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,D} {7,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-803.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,286,619,818,1246,1924,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.860571,'amu*angstrom^2'), symmetry=1, barrier=(19.7862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860595,'amu*angstrom^2'), symmetry=1, barrier=(19.7868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5331,'amu*angstrom^2'), symmetry=1, barrier=(35.2489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0557931,0.0951286,-0.000160821,1.39946e-07,-4.77786e-11,-96523.4,30.8142], Tmin=(100,'K'), Tmax=(800.001,'K')), NASAPolynomial(coeffs=[11.7634,0.0268665,-1.4597e-05,2.89929e-09,-2.03528e-13,-98085.5,-21.1111], Tmin=(800.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-803.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(COCFO) + radical(C=OCOJ)"""),
)

species(
    label = 'O=[C]C(=O)F(1543)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u1 p0 c0 {3,D} {4,S}
"""),
    E0 = (-325.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([233,496,705,1150,2014,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.06462,'amu*angstrom^2'), symmetry=1, barrier=(24.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0185,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90371,0.00648742,6.66758e-05,-1.78696e-07,1.33913e-10,-39126.2,10.0186], Tmin=(10,'K'), Tmax=(477.098,'K')), NASAPolynomial(coeffs=[4.98548,0.0142778,-1.08248e-05,3.66779e-09,-4.58866e-13,-39421.3,3.58924], Tmin=(477.098,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-325.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OD[C]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C[C](O)F(1553)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-416.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,280,501,1494,1531,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(2.02598,'amu*angstrom^2'), symmetry=1, barrier=(46.5814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0261,'amu*angstrom^2'), symmetry=1, barrier=(46.584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96882,0.00232875,9.96881e-05,-2.33543e-07,1.79738e-10,-50133.2,9.88029], Tmin=(10,'K'), Tmax=(332.028,'K')), NASAPolynomial(coeffs=[1.77552,0.0287515,-1.96799e-05,6.12833e-09,-7.19801e-13,-49987.5,18.0435], Tmin=(332.028,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-416.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC[C](O)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C=C([O])C(=O)C(=O)F(2436)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  O u1 p2 c0 {7,S}
3  O u0 p2 c0 {6,D}
4  O u1 p2 c0 {9,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,D} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u0 p0 c0 {4,S} {7,D} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-476.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,350,440,435,1725,286,619,818,1246,1924,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05437,'amu*angstrom^2'), symmetry=1, barrier=(24.2421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0571,'amu*angstrom^2'), symmetry=1, barrier=(24.3048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121386,0.0774108,-8.93685e-05,4.21166e-08,-5.60931e-12,-57128,28.0192], Tmin=(100,'K'), Tmax=(944.2,'K')), NASAPolynomial(coeffs=[22.2326,0.00243414,4.33263e-08,-4.405e-11,1.4745e-15,-62136.8,-81.7979], Tmin=(944.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-476.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsOsH) + group(COCFO) + radical(C=OCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(F)(C=O)C(=O)[C]=O(2375)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,D} {6,S} {9,S}
8  C u0 p0 c0 {4,D} {6,S} {10,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-380.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04027,'amu*angstrom^2'), symmetry=1, barrier=(23.9178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04034,'amu*angstrom^2'), symmetry=1, barrier=(23.9195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04024,'amu*angstrom^2'), symmetry=1, barrier=(23.9171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366868,0.0886179,-0.000155913,1.37812e-07,-4.6948e-11,-45604.9,29.5094], Tmin=(100,'K'), Tmax=(835.792,'K')), NASAPolynomial(coeffs=[10.9093,0.0234664,-1.26092e-05,2.46713e-09,-1.70659e-13,-46853.8,-16.3897], Tmin=(835.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(C(=O)F)C(=O)[C]=O(2437)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {3,D} {6,S} {9,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-413.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,486,617,768,1157,1926,1855,455,950,265.094,265.572,2771.93],'cm^-1')),
        HinderedRotor(inertia=(0.00242981,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712707,'amu*angstrom^2'), symmetry=1, barrier=(35.8996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158276,'amu*angstrom^2'), symmetry=1, barrier=(7.85193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.876417,0.0750328,-0.000122214,1.05111e-07,-3.5692e-11,-49602.5,31.7068], Tmin=(100,'K'), Tmax=(802.221,'K')), NASAPolynomial(coeffs=[9.77298,0.0234219,-1.21529e-05,2.38025e-09,-1.66138e-13,-50796.6,-7.79917], Tmin=(802.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCsFO) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)(C=O)[C](O)C(=O)F(2438)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {8,S} {12,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {7,S} {10,S}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {6,D} {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-716.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.165015,0.0991056,-0.000161875,1.35644e-07,-4.48198e-11,-86027.9,35.4458], Tmin=(100,'K'), Tmax=(781.34,'K')), NASAPolynomial(coeffs=[13.3389,0.0259684,-1.37795e-05,2.72387e-09,-1.91174e-13,-88015.9,-25.5887], Tmin=(781.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(COCsFO) + radical(C=OCOJ) + radical(C2CsJOH)"""),
)

species(
    label = '[O]C(F)=C([O])C(O)(F)C=O(2439)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {10,S}
3  O u0 p2 c0 {7,S} {12,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {7,S} {10,D}
9  C u0 p0 c0 {5,D} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {6,S} {8,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-786.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.263323,0.103516,-0.000176392,1.55722e-07,-5.3814e-11,-94420.9,31.9597], Tmin=(100,'K'), Tmax=(806.852,'K')), NASAPolynomial(coeffs=[11.6567,0.0310515,-1.68179e-05,3.3353e-09,-2.33826e-13,-95909.2,-20.2921], Tmin=(806.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-786.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(C(=O)F)C(O)(F)[C]=O(2440)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {12,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
8  C u0 p0 c0 {4,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {2,S} {5,D} {8,S}
10 C u1 p0 c0 {6,D} {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-733.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0955999,0.0910075,-0.000130579,9.36679e-08,-2.64077e-11,-88039.3,35.918], Tmin=(100,'K'), Tmax=(871.114,'K')), NASAPolynomial(coeffs=[15.2116,0.021597,-1.10584e-05,2.19757e-09,-1.56553e-13,-90672.8,-34.9223], Tmin=(871.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)CsOsH) + group(COCsFO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)([C]=O)C(O)C(=O)F(2441)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {12,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {5,D} {7,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-733.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0956011,0.0910075,-0.000130579,9.36679e-08,-2.64076e-11,-88039.3,35.918], Tmin=(100,'K'), Tmax=(871.118,'K')), NASAPolynomial(coeffs=[15.2116,0.021597,-1.10584e-05,2.19757e-09,-1.56553e-13,-90672.8,-34.9223], Tmin=(871.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=C(OF)C([O])C(=O)F(2442)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-404.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,1380,1390,370,380,2900,435,350,440,435,1725,486,617,768,1157,1926,3010,987.5,1337.5,450,1655,450.459,450.461,450.463,450.463,450.465],'cm^-1')),
        HinderedRotor(inertia=(0.000830772,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099409,'amu*angstrom^2'), symmetry=1, barrier=(14.314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099407,'amu*angstrom^2'), symmetry=1, barrier=(14.314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0201968,0.084752,-0.000102514,5.94635e-08,-1.33211e-11,-48525.9,36.042], Tmin=(100,'K'), Tmax=(1099.62,'K')), NASAPolynomial(coeffs=[19.3649,0.0143846,-6.52728e-06,1.27039e-09,-9.10372e-14,-52780.4,-59.1222], Tmin=(1099.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(COCsFO) + radical(C=OCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C([O])C(OF)C(=O)F(2443)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {3,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {4,S} {7,S} {10,D}
9  C u0 p0 c0 {1,S} {5,D} {7,S}
10 C u0 p0 c0 {6,S} {8,D} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-533.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,350,440,435,1725,486,617,768,1157,1926,3010,987.5,1337.5,450,1655,382.462,382.495,382.508,382.538,382.57],'cm^-1')),
        HinderedRotor(inertia=(0.177295,'amu*angstrom^2'), symmetry=1, barrier=(18.406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177293,'amu*angstrom^2'), symmetry=1, barrier=(18.4061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177324,'amu*angstrom^2'), symmetry=1, barrier=(18.4059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.638682,0.0894257,-0.000105258,5.71812e-08,-1.17048e-11,-64036.3,36.8314], Tmin=(100,'K'), Tmax=(1257.99,'K')), NASAPolynomial(coeffs=[25.0476,0.00488695,-1.04042e-06,1.41131e-10,-9.54328e-15,-70272.3,-92.0845], Tmin=(1257.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-533.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(COCsFO) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C([C]=O)C(F)(C=O)OF(2444)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-344.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,1855,455,950,256.726,256.826,256.921],'cm^-1')),
        HinderedRotor(inertia=(0.394483,'amu*angstrom^2'), symmetry=1, barrier=(18.4947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395056,'amu*angstrom^2'), symmetry=1, barrier=(18.4904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395272,'amu*angstrom^2'), symmetry=1, barrier=(18.4938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393798,'amu*angstrom^2'), symmetry=1, barrier=(18.4921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0561652,0.092321,-0.00013087,9.28189e-08,-2.59563e-11,-41312.4,35.3086], Tmin=(100,'K'), Tmax=(876.659,'K')), NASAPolynomial(coeffs=[15.2823,0.0228444,-1.19875e-05,2.40983e-09,-1.72939e-13,-43981.9,-36.144], Tmin=(876.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)(C=O)C([C]=O)OF(2445)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-344.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2782.5,750,1395,475,1775,1000,1855,455,950,256.726,256.826,256.921],'cm^-1')),
        HinderedRotor(inertia=(0.394483,'amu*angstrom^2'), symmetry=1, barrier=(18.4947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395056,'amu*angstrom^2'), symmetry=1, barrier=(18.4904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395272,'amu*angstrom^2'), symmetry=1, barrier=(18.4938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393798,'amu*angstrom^2'), symmetry=1, barrier=(18.4921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0561652,0.092321,-0.00013087,9.28189e-08,-2.59563e-11,-41312.4,35.3086], Tmin=(100,'K'), Tmax=(876.659,'K')), NASAPolynomial(coeffs=[15.2823,0.0228444,-1.19875e-05,2.40983e-09,-1.72939e-13,-43981.9,-36.144], Tmin=(876.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    E0 = (-215.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (111.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (181.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (227.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-80.4407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-73.1596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (147.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (153.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-207.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-151.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (41.6463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (87.8688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-56.9864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-24.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-113.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-93.1418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-93.3289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-93.3289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-29.3894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-46.0593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-191.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-136.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-105.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-135.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (12.1396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-73.0689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-4.43341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (41.7536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (26.5722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-61.4638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-140.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-140.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-87.9783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (127.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (74.1615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (168.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (193.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['O=CC(=O)F(335)', 'O=CC(=O)F(335)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C([O])C(=O)F(1893)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(88.9,'m^3/(mol*s)'), n=1.51, Ea=(338.188,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CO(13)', '[O]C(F)C([O])(F)C=O(1894)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(373.072,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(=COF)C([O])C(=O)F(2425)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(201.559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(O[CH]C(=O)F)C(=O)F(1907)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(134.931,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(F)(C=O)O[C](F)C=O(1911)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(142.212,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(6)', '[O]C=C(F)C([O])C(=O)F(2426)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(6)', '[O]C(F)([CH]C(=O)F)C=O(2427)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['O=CC1(F)OOC1C(=O)F(1916)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['O=CC(O)(F)C(=O)C(=O)F(2428)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(C(=O)F)C1(F)[CH]OO1(2429)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(257.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 255.4 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(F)(C=O)C1OO[C]1F(2430)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C1(F)[CH]OOC1C(=O)F(2305)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(158.385,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 155.9 to 158.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C1[C](F)OOC1(F)C=O(2254)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(190.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Endocyclic
Ea raised from 189.2 to 190.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C1OC1(F)C([O])C(=O)F(2431)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(F)(C=O)C1OC1([O])F(2432)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(122.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 122.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C1OC(C(=O)F)C1([O])F(2408)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C1C([O])(F)OC1(F)C=O(2324)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)(C=O)C(=O)[C](O)F(2433)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(258.339,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F(37)', '[O]C(C(=O)F)C(=O)C=O(2434)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(20.3257,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=CC(=O)F(335)', '[O][CH]C(=O)F(509)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(71.8791,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCO(15)', '[O]C(C(=O)F)C(=O)F(2235)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.04e+06,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(116.114,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CFO(51)', '[O]C(F)(C=O)C=O(2195)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.04e+06,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(83.4683,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(5)', '[O]C(F)(C=O)C(=O)C(=O)F(2435)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(22.3021,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][CH]C(=O)F(509)', '[O][CH]C(=O)F(509)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(7.10638,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['O=[C]C(=O)F(1543)', 'O=C[C](O)F(1553)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.01084e+12,'s^-1'), n=0.0883205, Ea=(142.303,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.021492990850914453, var=5.303057773824753, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-4R!H-R_7BrCClFIOPSSi->O',), comment="""Estimated from node Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_Ext-4R!H-R_7BrCClFIOPSSi->O"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', '[O]C=C([O])C(=O)C(=O)F(2436)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(318.903,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['HF(38)', '[O]C(F)(C=O)C(=O)[C]=O(2375)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(269.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HF(38)', '[O]C(C(=O)F)C(=O)[C]=O(2437)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(287.002,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(F)(C=O)[C](O)C(=O)F(2438)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(F)=C([O])C(O)(F)C=O(2439)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(C(=O)F)C(O)(F)[C]=O(2440)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    products = ['[O]C(F)([C]=O)C(O)C(=O)F(2441)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C=C(OF)C([O])C(=O)F(2442)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(98.2596,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=C([O])C(OF)C(=O)F(2443)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(174.081,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C([C]=O)C(F)(C=O)OF(2444)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(79.1809,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C(F)(C=O)C([C]=O)OF(2445)'],
    products = ['[O]C(C(=O)F)C([O])(F)C=O(1906)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #351',
    isomers = [
        '[O]C(C(=O)F)C([O])(F)C=O(1906)',
    ],
    reactants = [
        ('O=CC(=O)F(335)', 'O=CC(=O)F(335)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #351',
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

