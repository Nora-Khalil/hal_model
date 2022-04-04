species(
    label = '[O]C(C(=O)F)C(F)(F)[C]=CF(5257)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-663.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,136,307,446,511,682,757,1180,1185,486,617,768,1157,1926,615,860,1140,1343,3152,1685,370,238.573,238.573,238.574,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00296182,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15593,'amu*angstrom^2'), symmetry=1, barrier=(6.29798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683307,'amu*angstrom^2'), symmetry=1, barrier=(27.5986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0302886,0.0951656,-0.000136338,1.03265e-07,-3.1555e-11,-79650.6,37.3121], Tmin=(100,'K'), Tmax=(797.186,'K')), NASAPolynomial(coeffs=[12.6759,0.0317158,-1.69526e-05,3.42831e-09,-2.46437e-13,-81666.8,-20.8296], Tmin=(797.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-663.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = 'O=CC(=O)F(4234)',
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
    label = '[O]C(F)C(F)(F)[C]=CF(5325)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {3,S} {5,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-515.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,391,562,707,872,1109,1210,1289,3137,615,860,1140,1343,3152,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.280002,'amu*angstrom^2'), symmetry=1, barrier=(6.43779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280859,'amu*angstrom^2'), symmetry=1, barrier=(6.45751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646887,0.0801067,-0.000120732,9.74801e-08,-3.16857e-11,-61913.5,28.5357], Tmin=(100,'K'), Tmax=(751.816,'K')), NASAPolynomial(coeffs=[10.7127,0.0265596,-1.39112e-05,2.77136e-09,-1.96795e-13,-63427.3,-17.156], Tmin=(751.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFH) + radical(O2sj(Cs-CsF1sH)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = '[O]C(C(=O)F)C([CH]F)=C(F)F(5249)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {7,S} {10,S} {11,D}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 C u1 p0 c0 {2,S} {8,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-742.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269151,0.0881154,-0.000112605,7.54746e-08,-2.04419e-11,-89137.5,36.7913], Tmin=(100,'K'), Tmax=(894.936,'K')), NASAPolynomial(coeffs=[13.2379,0.0301492,-1.54461e-05,3.09615e-09,-2.22575e-13,-91458.7,-24.3354], Tmin=(894.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHH) + group(Cds-CdsCsCs) + group(COCsFO) + group(CdCFF) + radical(C=OCOJ) + radical(Csj(Cd-CsCd)(F1s)(H))"""),
)

species(
    label = 'O=C[C](F)OC(F)(F)[C]=CF(5260)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
8  C u1 p0 c0 {3,S} {5,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {12,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-715.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([140,424,499,621,667,843,876,1082,280,501,1494,1531,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,180,180,180,1590.28],'cm^-1')),
        HinderedRotor(inertia=(0.308241,'amu*angstrom^2'), symmetry=1, barrier=(7.08708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308375,'amu*angstrom^2'), symmetry=1, barrier=(7.09014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26114,'amu*angstrom^2'), symmetry=1, barrier=(51.988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26158,'amu*angstrom^2'), symmetry=1, barrier=(51.9983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3688.53,'J/mol'), sigma=(5.96387,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=576.14 K, Pc=39.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.388823,0.108847,-0.00018847,1.73094e-07,-6.18421e-11,-85900.5,33.5812], Tmin=(100,'K'), Tmax=(819.238,'K')), NASAPolynomial(coeffs=[9.31038,0.039487,-2.11876e-05,4.18195e-09,-2.92192e-13,-86751.3,-6.77162], Tmin=(819.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(CsCOF1sO2s) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
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
    label = 'O=C(F)[CH]C(F)(F)[C]=CF(5348)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u1 p0 c0 {6,S} {8,S} {11,S}
8  C u0 p0 c0 {3,S} {5,D} {7,S}
9  C u0 p0 c0 {4,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-575.105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,3025,407.5,1350,352.5,611,648,830,1210,1753,615,860,1140,1343,3152,1685,370,189.575,190.647,2456.43],'cm^-1')),
        HinderedRotor(inertia=(0.272591,'amu*angstrom^2'), symmetry=1, barrier=(6.9964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31945,'amu*angstrom^2'), symmetry=1, barrier=(33.8249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46064,'amu*angstrom^2'), symmetry=1, barrier=(33.8156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.359045,0.0877677,-0.000126492,9.81967e-08,-3.10049e-11,-69045,30.9214], Tmin=(100,'K'), Tmax=(770.114,'K')), NASAPolynomial(coeffs=[11.1745,0.0315883,-1.70615e-05,3.45986e-09,-2.48954e-13,-70710.8,-18.4313], Tmin=(770.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-575.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(CCJC=O) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
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
    label = '[O]C([C](F)F)C(=O)F(4487)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {1,S} {2,S} {6,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-591.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,190,488,555,1236,1407,486,617,768,1157,1926,378.555,378.585,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0876046,'amu*angstrom^2'), symmetry=1, barrier=(8.91333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259265,'amu*angstrom^2'), symmetry=1, barrier=(26.381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43774,0.0608417,-8.50182e-05,6.14589e-08,-1.77788e-11,-71047.5,27.7795], Tmin=(100,'K'), Tmax=(843.24,'K')), NASAPolynomial(coeffs=[10.3543,0.0185442,-9.77556e-06,1.97086e-09,-1.41659e-13,-72551.2,-13.7173], Tmin=(843.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-591.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'O=C(F)C1OC(=CF)C1(F)F(5267)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {8,S} {11,D}
10 C u0 p0 c0 {3,S} {6,D} {7,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1045.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.140035,0.0758762,-4.96867e-05,-7.50456e-09,1.3059e-11,-125550,28.405], Tmin=(100,'K'), Tmax=(945.047,'K')), NASAPolynomial(coeffs=[24.609,0.00835338,-1.60499e-06,2.62536e-10,-2.3074e-14,-131890,-98.3919], Tmin=(945.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1045.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsOs) + group(COCsFO) + group(CdCFH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'O=C(F)C(=O)C(F)(F)C=CF(5349)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  C u0 p0 c0 {7,S} {11,D} {12,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-1065.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0233108,0.0982609,-0.000140671,9.61765e-08,-2.0423e-11,-127960,31.0377], Tmin=(100,'K'), Tmax=(596.94,'K')), NASAPolynomial(coeffs=[12.0312,0.0336239,-1.8018e-05,3.61092e-09,-2.56848e-13,-129676,-23.0602], Tmin=(596.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1065.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(COCFO) + group(CdCFH)"""),
)

species(
    label = 'FC=[C]C(F)(F)C1OO[C]1F(5350)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u1 p0 c0 {3,S} {6,S} {7,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-360.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.132704,0.0995716,-0.000158168,1.37029e-07,-4.77782e-11,-43192.6,34.757], Tmin=(100,'K'), Tmax=(756.821,'K')), NASAPolynomial(coeffs=[10.8687,0.0353674,-1.89087e-05,3.78085e-09,-2.682e-13,-44684.3,-14.1068], Tmin=(756.821,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(CsCCFF) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = '[O]C1[C](F)OC(=CF)C1(F)F(5351)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {5,S} {7,S} {11,D}
10 C u1 p0 c0 {3,S} {5,S} {8,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-742.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.31389,0.0820167,-7.19411e-05,2.22174e-08,7.52783e-14,-89122.3,30.3664], Tmin=(100,'K'), Tmax=(1033.18,'K')), NASAPolynomial(coeffs=[23.5064,0.0130269,-5.50761e-06,1.11384e-09,-8.43171e-14,-95284.4,-91.3307], Tmin=(1033.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFF) + group(CsCFHO) + group(Cds-CdsCsOs) + group(CdCFH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)OC1C(F)(F)[C]=CF(5352)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
10 C u0 p0 c0 {4,S} {11,D} {13,S}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-541.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.240248,0.0905853,-0.000121415,8.80988e-08,-2.61867e-11,-64957.2,31.7575], Tmin=(100,'K'), Tmax=(813.749,'K')), NASAPolynomial(coeffs=[11.6067,0.0347157,-1.84332e-05,3.73451e-09,-2.69407e-13,-66807.1,-20.7369], Tmin=(813.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-541.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFOO) + group(CsCCFF) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = '[O]C1C([O])(F)C(=CF)C1(F)F(5353)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
10 C u0 p0 c0 {8,S} {9,S} {11,D}
11 C u0 p0 c0 {4,S} {10,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-632.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117869,0.0935693,-0.000118376,7.80947e-08,-2.1006e-11,-75946.5,26.2213], Tmin=(100,'K'), Tmax=(895.78,'K')), NASAPolynomial(coeffs=[13.4612,0.0339843,-1.85961e-05,3.83292e-09,-2.79836e-13,-78337,-36.6833], Tmin=(895.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-632.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCCFF) + group(Cds-CdsCsCs) + group(CdCFH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=C([C](O)F)C(F)(F)[C]=CF(5354)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {6,D} {7,S} {9,S}
9  C u1 p0 c0 {3,S} {5,S} {8,S}
10 C u0 p0 c0 {4,S} {11,D} {12,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-745.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,280,501,1494,1531,615,860,1140,1343,3152,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.58687,'amu*angstrom^2'), symmetry=1, barrier=(36.4852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58678,'amu*angstrom^2'), symmetry=1, barrier=(36.4831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960983,'amu*angstrom^2'), symmetry=1, barrier=(22.0949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961651,'amu*angstrom^2'), symmetry=1, barrier=(22.1102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.608022,0.115599,-0.000209599,1.93489e-07,-6.81348e-11,-89480.2,34.1756], Tmin=(100,'K'), Tmax=(846.723,'K')), NASAPolynomial(coeffs=[9.84683,0.0373576,-1.98785e-05,3.8701e-09,-2.66619e-13,-90216.4,-8.4161], Tmin=(846.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-745.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(CdCFH) + radical(CsCOF1sO2s) + radical(Cds_S)"""),
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
    label = 'O=CC(F)(F)[C]=CF(3294)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {4,D} {5,S} {9,S}
7  C u0 p0 c0 {3,S} {8,D} {10,S}
8  C u1 p0 c0 {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-455.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,615,860,1140,1343,3152,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.26344,'amu*angstrom^2'), symmetry=1, barrier=(29.0491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26217,'amu*angstrom^2'), symmetry=1, barrier=(29.0199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07634,0.0717469,-0.000119864,1.08093e-07,-3.82637e-11,-54709.8,24.2324], Tmin=(100,'K'), Tmax=(812.992,'K')), NASAPolynomial(coeffs=[7.79931,0.0264797,-1.38546e-05,2.7214e-09,-1.90074e-13,-55400.1,-4.3327], Tmin=(812.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-455.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(Cds_S)"""),
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
    label = 'O=C(F)C(=O)C(F)(F)[C]=CF(5355)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u0 p0 c0 {5,D} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 C u0 p0 c0 {4,S} {11,D} {12,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-827.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,286,619,818,1246,1924,615,860,1140,1343,3152,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15539,'amu*angstrom^2'), symmetry=1, barrier=(26.5647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1548,'amu*angstrom^2'), symmetry=1, barrier=(26.5512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15563,'amu*angstrom^2'), symmetry=1, barrier=(26.5702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425309,0.108819,-0.000192678,1.73835e-07,-6.06836e-11,-99339.3,33.0614], Tmin=(100,'K'), Tmax=(822.02,'K')), NASAPolynomial(coeffs=[11.3964,0.032105,-1.76748e-05,3.5073e-09,-2.45133e-13,-100634,-17.7109], Tmin=(822.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-827.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(COCFO) + group(CdCFH) + radical(Cds_S)"""),
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
    label = '[O]C(C(=O)F)C(F)=C=CF(5356)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {10,D}
8  C u0 p0 c0 {2,S} {5,D} {6,S}
9  C u0 p0 c0 {3,S} {10,D} {12,S}
10 C u0 p0 c0 {7,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-504.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,145,326,398,834,1303,486,617,768,1157,1926,113,247,382,1207,3490,540,610,2055,873.164,875.045,876.694,2871.62],'cm^-1')),
        HinderedRotor(inertia=(0.0718076,'amu*angstrom^2'), symmetry=1, barrier=(4.39584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705624,'amu*angstrom^2'), symmetry=1, barrier=(4.39846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586258,0.0812342,-0.000113217,8.46996e-08,-2.56177e-11,-60588.1,32.952], Tmin=(100,'K'), Tmax=(805.094,'K')), NASAPolynomial(coeffs=[11.2629,0.0281881,-1.43833e-05,2.85786e-09,-2.03551e-13,-62307.2,-16.2416], Tmin=(805.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-504.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CdCddCF) + group(COCsFO) + group(CdCddFH) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(398)',
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
    label = 'C#CC(F)(F)C([O])C(=O)F(5357)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {5,D} {6,S}
9  C u0 p0 c0 {7,S} {10,T}
10 C u0 p0 c0 {9,T} {12,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-564.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,154,355,414,641,686,1150,1196,486,617,768,1157,1926,2175,525,750,770,3400,2100,180,229.569,2555.64],'cm^-1')),
        HinderedRotor(inertia=(0.203636,'amu*angstrom^2'), symmetry=1, barrier=(15.3525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199055,'amu*angstrom^2'), symmetry=1, barrier=(15.3533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6497,'amu*angstrom^2'), symmetry=1, barrier=(49.024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.228522,0.0874893,-0.000122557,8.79336e-08,-2.49264e-11,-67806.8,31.8097], Tmin=(100,'K'), Tmax=(865.449,'K')), NASAPolynomial(coeffs=[14.1256,0.0232574,-1.1228e-05,2.17332e-09,-1.52514e-13,-70212.2,-33.227], Tmin=(865.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-564.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(COCsFO) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)C#CF(5358)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u0 p0 c0 {8,S} {11,T}
11 C u0 p0 c0 {4,S} {10,T}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-664.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,154,355,414,641,686,1150,1196,486,617,768,1157,1926,2175,525,239,401,1367,233.873,233.908,1099.99,1099.99],'cm^-1')),
        HinderedRotor(inertia=(0.188389,'amu*angstrom^2'), symmetry=1, barrier=(7.31158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188366,'amu*angstrom^2'), symmetry=1, barrier=(7.31145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00308172,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.117209,0.0981222,-0.000158553,1.33431e-07,-4.44267e-11,-79796.6,35.3617], Tmin=(100,'K'), Tmax=(784.986,'K')), NASAPolynomial(coeffs=[12.6561,0.0279139,-1.4611e-05,2.87568e-09,-2.01477e-13,-81644.2,-22.165], Tmin=(784.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-664.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(COCsFO) + group(Ct-CtCs) + group(CtCF) + radical(C=OCOJ)"""),
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
    label = 'O=C(F)C(=O)C(F)=[C][CH]F(5359)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,D} {7,S} {8,S}
7  C u0 p0 c0 {1,S} {6,S} {10,D}
8  C u0 p0 c0 {2,S} {5,D} {6,S}
9  C u1 p0 c0 {3,S} {10,S} {11,S}
10 C u1 p0 c0 {7,D} {9,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-485.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,244,384,691,1241,286,619,818,1246,1924,234,589,736,816,1240,3237,1685,370,277.707,279.44,1497.85],'cm^-1')),
        HinderedRotor(inertia=(0.112192,'amu*angstrom^2'), symmetry=1, barrier=(6.12223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304876,'amu*angstrom^2'), symmetry=1, barrier=(48.5366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.887694,'amu*angstrom^2'), symmetry=1, barrier=(48.5338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.616532,0.0811913,-0.000119159,9.46399e-08,-3.04998e-11,-58318.7,29.5606], Tmin=(100,'K'), Tmax=(755.983,'K')), NASAPolynomial(coeffs=[10.506,0.0288641,-1.53308e-05,3.07702e-09,-2.19785e-13,-59813.9,-15.3836], Tmin=(755.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-485.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]C(=O)C(F)(F)[C]=CF(5337)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {10,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {4,D} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {9,D} {11,S}
9  C u1 p0 c0 {6,S} {8,D}
10 C u1 p0 c0 {5,D} {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-403.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,615,860,1140,1343,3152,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16242,'amu*angstrom^2'), symmetry=1, barrier=(26.7264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16411,'amu*angstrom^2'), symmetry=1, barrier=(26.7653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16288,'amu*angstrom^2'), symmetry=1, barrier=(26.7369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.102227,0.102157,-0.0001872,1.70915e-07,-5.95064e-11,-48421.2,31.7142], Tmin=(100,'K'), Tmax=(844.589,'K')), NASAPolynomial(coeffs=[10.4973,0.0287853,-1.5735e-05,3.08679e-09,-2.13249e-13,-49385.2,-12.7392], Tmin=(844.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(CdCFH) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C(F)[C](O)C(F)(F)[C]=CF(5360)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u1 p0 c0 {5,S} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 C u0 p0 c0 {4,S} {11,D} {12,S}
11 C u1 p0 c0 {7,S} {10,D}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-730.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.572948,0.109591,-0.000180878,1.54658e-07,-5.23582e-11,-87701.3,38.1015], Tmin=(100,'K'), Tmax=(771.059,'K')), NASAPolynomial(coeffs=[13.5119,0.0308975,-1.68452e-05,3.37082e-09,-2.3837e-13,-89706.1,-25.1032], Tmin=(771.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-730.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C2CsJOH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH))"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)C=[C]F(5361)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 C u1 p0 c0 {4,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-672.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.091246,0.0976849,-0.000143783,1.11044e-07,-3.43845e-11,-80742.4,36.6904], Tmin=(100,'K'), Tmax=(789.058,'K')), NASAPolynomial(coeffs=[13.1457,0.0305801,-1.62119e-05,3.25754e-09,-2.32932e-13,-82831.3,-24.0339], Tmin=(789.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-672.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(C=OCOJ) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = '[O]C(F)=C([O])C(F)(F)C=CF(5362)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {9,S}
6  O u1 p2 c0 {11,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u0 p0 c0 {5,S} {7,S} {11,D}
10 C u0 p0 c0 {3,S} {8,D} {13,S}
11 C u0 p0 c0 {4,S} {6,S} {9,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-814.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.585661,0.115462,-0.000209367,1.96282e-07,-7.06143e-11,-97825,33.0482], Tmin=(100,'K'), Tmax=(829.903,'K')), NASAPolynomial(coeffs=[9.04627,0.0404558,-2.21391e-05,4.38244e-09,-3.05714e-13,-98439.4,-5.69446], Tmin=(829.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(CdCFH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C(F)C(O)C(F)(F)[C]=[C]F(5363)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {13,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {4,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-656.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,136,307,446,511,682,757,1180,1185,486,617,768,1157,1926,1685,370,167,640,1190,210.624,222.328],'cm^-1')),
        HinderedRotor(inertia=(0.435843,'amu*angstrom^2'), symmetry=1, barrier=(12.5185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116001,'amu*angstrom^2'), symmetry=1, barrier=(5.08468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00439013,'amu*angstrom^2'), symmetry=1, barrier=(0.141851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0302309,'amu*angstrom^2'), symmetry=1, barrier=(12.5914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.218654,0.100836,-0.000149247,1.14086e-07,-3.48353e-11,-78873.9,37.9152], Tmin=(100,'K'), Tmax=(800.723,'K')), NASAPolynomial(coeffs=[13.9984,0.0298066,-1.61737e-05,3.27995e-09,-2.35779e-13,-81150.4,-27.5127], Tmin=(800.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'O=C(F)C(OF)C(F)=[C][CH]F(5364)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {11,D}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 C u1 p0 c0 {3,S} {11,S} {13,S}
11 C u1 p0 c0 {8,D} {10,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-428.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,246,474,533,1155,486,617,768,1157,1926,234,589,736,816,1240,3237,1685,370,314.884,314.896,314.925,1747.19],'cm^-1')),
        HinderedRotor(inertia=(0.00170207,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143619,'amu*angstrom^2'), symmetry=1, barrier=(10.1114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411308,'amu*angstrom^2'), symmetry=1, barrier=(28.9539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411681,'amu*angstrom^2'), symmetry=1, barrier=(28.9535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.133176,0.0984201,-0.000134576,9.42267e-08,-2.64328e-11,-51446,35.9509], Tmin=(100,'K'), Tmax=(867.68,'K')), NASAPolynomial(coeffs=[14.7127,0.0299808,-1.62618e-05,3.32173e-09,-2.40784e-13,-54022.3,-33.5646], Tmin=(867.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-428.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHH) + group(CdCsCdF) + group(Cds-CdsCsH) + group(COCsFO) + radical(Csj(Cd-CdH)(F1s)(H)) + radical(Cdj(Cs-F1sHH)(Cd-CsF1s))"""),
)

species(
    label = '[O]C([C](F)C(F)=CF)C(=O)F(5365)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {12,S}
8  C u1 p0 c0 {1,S} {7,S} {9,S}
9  C u0 p0 c0 {2,S} {8,S} {11,D}
10 C u0 p0 c0 {3,S} {6,D} {7,S}
11 C u0 p0 c0 {4,S} {9,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-730.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0418831,0.0933784,-0.000124659,8.61215e-08,-2.38467e-11,-87718.4,36.1764], Tmin=(100,'K'), Tmax=(879.189,'K')), NASAPolynomial(coeffs=[14.1577,0.0291537,-1.50795e-05,3.02664e-09,-2.17454e-13,-90200.4,-30.1063], Tmin=(879.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-730.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(COCsFO) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(C=OCOJ) + radical(CsCdCsF1s)"""),
)

species(
    label = 'O=[C]C(OF)C(F)(F)[C]=CF(5366)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {10,D} {13,S}
10 C u1 p0 c0 {8,S} {9,D}
11 C u1 p0 c0 {6,D} {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-358.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,136,307,446,511,682,757,1180,1185,615,860,1140,1343,3152,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00520297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00680427,'amu*angstrom^2'), symmetry=1, barrier=(0.156444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.686555,'amu*angstrom^2'), symmetry=1, barrier=(15.7853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68534,'amu*angstrom^2'), symmetry=1, barrier=(15.7573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.299837,0.102099,-0.000146897,1.07131e-07,-3.10427e-11,-42987.9,37.7839], Tmin=(100,'K'), Tmax=(843.939,'K')), NASAPolynomial(coeffs=[15.335,0.0279955,-1.51887e-05,3.09013e-09,-2.22987e-13,-45626.9,-34.9925], Tmin=(843.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CdCFH) + radical(Cdj(Cs-F1sF1sCs)(Cd-F1sH)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C([C]=O)C(F)(F)C(F)=CF(5367)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {11,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {13,S}
11 C u1 p0 c0 {6,D} {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-675.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.438949,0.104871,-0.000163433,1.30486e-07,-4.114e-11,-81103.2,37.3787], Tmin=(100,'K'), Tmax=(779.949,'K')), NASAPolynomial(coeffs=[14.7059,0.0271905,-1.40194e-05,2.75754e-09,-1.93755e-13,-83465.4,-31.9209], Tmin=(779.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-675.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C]C(F)(F)C(OF)C(=O)F(5368)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 C u1 p0 c0 {8,S} {11,D}
11 C u1 p0 c0 {10,D} {13,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-350.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,136,307,446,511,682,757,1180,1185,486,617,768,1157,1926,1685,370,3120,650,792.5,1650,180,756.653],'cm^-1')),
        HinderedRotor(inertia=(0.272809,'amu*angstrom^2'), symmetry=1, barrier=(6.27241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845227,'amu*angstrom^2'), symmetry=1, barrier=(19.4334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846473,'amu*angstrom^2'), symmetry=1, barrier=(19.4621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582111,'amu*angstrom^2'), symmetry=1, barrier=(13.3839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.312919,0.10298,-0.000147539,1.07354e-07,-3.11421e-11,-42004,37.1998], Tmin=(100,'K'), Tmax=(841.615,'K')), NASAPolynomial(coeffs=[15.1911,0.0292935,-1.62109e-05,3.3258e-09,-2.41177e-13,-44613.7,-34.9247], Tmin=(841.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-350.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + group(Cds-CdsCsH) + group(COCsFO) + group(Cds-CdsHH) + radical(Cdj(Cs-F1sF1sCs)(Cd-HH)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(F)C(F)(F)C([O])C(=O)F(5369)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u0 p0 c0 {4,S} {6,D} {8,S}
11 C u1 p0 c0 {9,D} {13,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {11,S}
"""),
    E0 = (-664.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,1380,1390,370,380,2900,435,246,474,533,1155,486,617,768,1157,1926,3120,650,792.5,1650,304.654,304.659,304.664,4000],'cm^-1')),
        HinderedRotor(inertia=(0.131973,'amu*angstrom^2'), symmetry=1, barrier=(8.69239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131972,'amu*angstrom^2'), symmetry=1, barrier=(8.69239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0222493,'amu*angstrom^2'), symmetry=1, barrier=(29.6917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (170.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20951,0.0998388,-0.000146184,1.1067e-07,-3.33748e-11,-79766.9,36.5234], Tmin=(100,'K'), Tmax=(811.598,'K')), NASAPolynomial(coeffs=[14.1755,0.0289393,-1.51427e-05,3.02606e-09,-2.15764e-13,-82101.8,-29.8727], Tmin=(811.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-664.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(COCsFO) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cdj(Cd-CsF1s)(H))"""),
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
    E0 = (-232.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (167.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-106.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-95.6213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (98.9744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (244.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-224.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-143.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (70.9038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-181.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-110.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-151.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-50.7933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-212.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-127.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-158.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (49.3594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-133.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (6.41887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-21.7895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (20.3264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-54.8278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (19.9329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-78.4288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-28.1995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-90.3599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-192.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (129.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-81.8281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (176.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (4.43766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (157.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-50.9193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', '[O]C(F)C(F)(F)[C]=CF(5325)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(370.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C(C(=O)F)C([CH]F)=C(F)F(5249)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.07321e+10,'s^-1'), n=0.899, Ea=(125.897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-R!HR!H)CJ;CJ;C] for rate rule [cCs(-R!HR!H)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['O=C[C](F)OC(F)(F)[C]=CF(5260)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(136.716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C(F)[CH]C(F)(F)[C]=CF(5348)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=CF-2(1206)', '[O]C([C](F)F)C(=O)F(4487)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['O=C(F)C1OC(=CF)C1(F)F(5267)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['O=C(F)C(=O)C(F)(F)C=CF(5349)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['FC=[C]C(F)(F)C1OO[C]1F(5350)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C1[C](F)OC(=CF)C1(F)F(5351)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.67521e+10,'s^-1'), n=0.355, Ea=(50.9402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cddouble] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_cddouble]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C1(F)OC1C(F)(F)[C]=CF(5352)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(122.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 122.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C1C([O])(F)C(=CF)C1(F)F(5353)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_cddouble]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C([C](O)F)C(F)(F)[C]=CF(5354)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(263.408,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=CC(=O)F(4234)', 'F[CH][C]=C(F)F(2842)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.36858e-17,'m^3/(mol*s)'), n=6.25044, Ea=(39.9792,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.6646553062382831, var=1.4576930825215244, Tref=1000.0, N=48, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing_Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_Ext-3R-R_Ext-1R!H-R_N-8R!H-inRing_Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CFO(51)', 'O=CC(F)(F)[C]=CF(3294)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(87.7051,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'O=C(F)C(=O)C(F)(F)[C]=CF(5355)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(25.3822,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', '[O]C(C(=O)F)C(F)=C=CF(5356)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(50.1565,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][CH]C(=O)F(398)', 'FC=C=C(F)F(1325)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(14.6253,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F(37)', 'C#CC(F)(F)C([O])C(=O)F(5357)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(67.356,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(C(=O)F)C(F)(F)C#CF(5358)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(10.8869,'m^3/(mol*s)'), n=2.06837, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.012984286287768347, var=0.42775789296846467, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][CH]C(=O)F(398)', 'F[CH][C]=C(F)F(2842)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(4.42411,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'O=C(F)C(=O)C(F)=[C][CH]F(5359)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(281.091,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'O=[C]C(=O)C(F)(F)[C]=CF(5337)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(273.731,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['O=C(F)[C](O)C(F)(F)[C]=CF(5360)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C(C(=O)F)C(F)(F)C=[C]F(5361)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C(F)=C([O])C(F)(F)C=CF(5362)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C(F)C(O)C(F)(F)[C]=[C]F(5363)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;XH_out] for rate rule [R5HJ_1;Cd_rad_out_single;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=C(F)C(OF)C(F)=[C][CH]F(5364)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(127.249,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C([C](F)C(F)=CF)C(=O)F(5365)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(150.509,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=[C]C(OF)C(F)(F)[C]=CF(5366)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    products = ['[O]C([C]=O)C(F)(F)C(F)=CF(5367)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(236.775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=[C]C(F)(F)C(OF)C(=O)F(5368)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(76.8185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(F)C(F)(F)C([O])C(=O)F(5369)'],
    products = ['[O]C(C(=O)F)C(F)(F)[C]=CF(5257)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(182.461,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #1536',
    isomers = [
        '[O]C(C(=O)F)C(F)(F)[C]=CF(5257)',
    ],
    reactants = [
        ('O=CC(=O)F(4234)', 'FC=C=C(F)F(1325)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1536',
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

