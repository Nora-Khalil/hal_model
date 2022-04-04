species(
    label = 'F[C]=CC1(F)[CH]O1(2817)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {8,S}
6 C u0 p0 c0 {4,S} {7,D} {9,S}
7 C u1 p0 c0 {2,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (5.69222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,376.683,377.053,377.179,377.675],'cm^-1')),
        HinderedRotor(inertia=(0.214823,'amu*angstrom^2'), symmetry=1, barrier=(21.6684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61816,0.0521954,-5.80947e-05,3.21258e-08,-7.02052e-12,770.605,21.228], Tmin=(100,'K'), Tmax=(1111.73,'K')), NASAPolynomial(coeffs=[12.1961,0.0141356,-6.74237e-06,1.33133e-09,-9.5586e-14,-1581.35,-30.9249], Tmin=(1111.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.69222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s) + radical(CCsJO) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'C2HF(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (95.331,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0062,'amu')),
        LinearRotor(inertia=(51.6236,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([429.793,429.793,596.357,596.357,1107.96,2365.05,3506.88],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1870.76,'J/mol'), sigma=(4.25,'angstroms'), dipoleMoment=(1,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4498,0.0030263,3.99146e-05,-8.9615e-08,5.74336e-11,11468.6,5.90915], Tmin=(10,'K'), Tmax=(555.749,'K')), NASAPolynomial(coeffs=[4.23833,0.0086714,-5.87678e-06,1.96876e-09,-2.53031e-13,11206.2,0.995297], Tmin=(555.749,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(95.331,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C#CF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CO1(1231)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-29.4903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,326,540,652,719,1357,788.854,788.882],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3010.76,'J/mol'), sigma=(4.86933,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=470.27 K, Pc=59.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.14414,0.0129866,6.60907e-06,-2.08037e-08,9.56931e-12,-3510.75,9.45209], Tmin=(100,'K'), Tmax=(957.696,'K')), NASAPolynomial(coeffs=[8.89006,0.00340757,-9.7292e-07,1.96224e-10,-1.667e-14,-5272.61,-21.4727], Tmin=(957.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.4903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + ring(Cyclopropene)"""),
)

species(
    label = 'F[C]=CC1O[C]1F(2818)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5 C u1 p0 c0 {1,S} {3,S} {4,S}
6 C u0 p0 c0 {4,S} {7,D} {9,S}
7 C u1 p0 c0 {2,S} {6,D}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (35.4612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,798.092,1149.53,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00379227,'amu*angstrom^2'), symmetry=1, barrier=(1.75346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3271.44,'J/mol'), sigma=(5.28954,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=510.99 K, Pc=50.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96402,0.0432867,-4.10789e-05,1.96399e-08,-3.75993e-12,4339.5,23.9355], Tmin=(100,'K'), Tmax=(1250.83,'K')), NASAPolynomial(coeffs=[10.9457,0.0145644,-6.63484e-06,1.28192e-09,-9.07465e-14,2092.6,-21.4058], Tmin=(1250.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.4612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'FC1=CC2(F)OC12(3086)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
6 C u0 p0 c0 {2,S} {4,S} {7,D}
7 C u0 p0 c0 {5,S} {6,D} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-250.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76628,0.0397831,-1.86277e-05,-6.18196e-09,4.91697e-12,-30036.3,14.5527], Tmin=(100,'K'), Tmax=(1112.21,'K')), NASAPolynomial(coeffs=[14.3569,0.0127358,-6.74169e-06,1.43386e-09,-1.08171e-13,-33964.8,-52.599], Tmin=(1112.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFO) + group(CdCsCdF) + group(Cds-CdsCsH) + polycyclic(s2_3_4_ene_1)"""),
)

species(
    label = 'O=CC(F)=C[C]F(3094)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {4,D} {7,S} {8,S}
6 C u0 p0 c0 {3,D} {4,S} {9,S}
7 C u2 p0 c0 {2,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-119.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([288,410,724,839,1320,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,180,180,1391.07],'cm^-1')),
        HinderedRotor(inertia=(0.113488,'amu*angstrom^2'), symmetry=1, barrier=(2.60931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03986,'amu*angstrom^2'), symmetry=1, barrier=(23.9084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17505,0.045907,-5.11787e-05,3.341e-08,-9.67677e-12,-14331.1,18.7068], Tmin=(100,'K'), Tmax=(803.475,'K')), NASAPolynomial(coeffs=[6.11536,0.0262905,-1.45565e-05,3.02326e-09,-2.2189e-13,-14964.2,0.559316], Tmin=(803.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFHH) + group(Cds-CdsCsH) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)H) + radical(CsCCl_triplet)"""),
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
    label = 'F[C]=CC1=CO1(3095)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {8,S}
5 C u0 p0 c0 {3,S} {6,D} {7,S}
6 C u1 p0 c0 {1,S} {5,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (266.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3010,987.5,1337.5,450,1655,167,640,1190,367.699,369.672,369.898,370.004,370.275,370.5,371.098],'cm^-1')),
        HinderedRotor(inertia=(0.00124301,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0563,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37809,0.0515362,-6.20075e-05,3.52807e-08,-7.57414e-12,32129,17.1825], Tmin=(100,'K'), Tmax=(1223.44,'K')), NASAPolynomial(coeffs=[14.9962,0.00427507,-7.06957e-07,4.86503e-11,-1.06975e-15,29001.6,-50.4261], Tmin=(1223.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(CdCFH) + ring(Cyclopropene) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = '[CH]=[C]F(2014)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u1 p0 c0 {3,D} {4,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (350.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([167,640,1190,1142.58,1502.03,3807.5],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.55591,0.0076799,-5.30098e-07,-4.54651e-09,2.16658e-12,42166.8,9.46363], Tmin=(100,'K'), Tmax=(1043.42,'K')), NASAPolynomial(coeffs=[5.86818,0.00351556,-1.29996e-06,2.62251e-10,-1.98915e-14,41428.4,-3.016], Tmin=(1043.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(CdCdF1s)"""),
)

species(
    label = 'F[C]1[CH]O1(2828)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u1 p0 c0 {2,S} {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (97.7552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,661.567,661.569,2935.1],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18921,0.0134531,-5.29464e-07,-1.06374e-08,5.3494e-12,11790.3,13.5601], Tmin=(100,'K'), Tmax=(993.343,'K')), NASAPolynomial(coeffs=[8.19763,0.00359113,-1.19992e-06,2.57126e-10,-2.11232e-14,10286.8,-13.1285], Tmin=(993.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.7552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(CsCsF1sO2s)"""),
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
    label = 'FC#CC1(F)[CH]O1(3096)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {8,S}
6 C u0 p0 c0 {4,S} {7,T}
7 C u0 p0 c0 {2,S} {6,T}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (13.5497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,2175,525,239,401,1367,249.072,249.251,249.574,249.762],'cm^-1')),
        HinderedRotor(inertia=(0.00270216,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.047,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76677,0.0503091,-6.34567e-05,4.03918e-08,-1.01288e-11,1709.11,19.2903], Tmin=(100,'K'), Tmax=(976.25,'K')), NASAPolynomial(coeffs=[10.893,0.0129164,-6.00353e-06,1.15805e-09,-8.17471e-14,-72.7789,-24.5189], Tmin=(976.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.5497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(CtCF) + ring(Cs-Cs(F)-O2s) + radical(CCsJO)"""),
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
    label = 'F[C]=C=C1[CH]O1(3097)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,S} {5,D}
4 C u1 p0 c0 {2,S} {3,S} {7,S}
5 C u0 p0 c0 {3,D} {6,D}
6 C u1 p0 c0 {1,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (316.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,540,610,2055,137,207,812,664.717,664.907,665.019,665.094,665.171,665.186,665.223],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0483,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22912,0.0211869,3.5459e-05,-7.77761e-08,3.5923e-11,38104.6,15.4478], Tmin=(100,'K'), Tmax=(918.033,'K')), NASAPolynomial(coeffs=[19.2793,-0.00678804,5.49255e-06,-1.05992e-09,6.62142e-14,33022.4,-75.981], Tmin=(918.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(CdCddFH) + group(Cdd-CdsCds) + ring(methyleneoxirane) + radical(C=CCJO) + radical(Cdj(Cdd-Cd)(F1s))"""),
)

species(
    label = 'FC=[C]C1(F)[CH]O1(3098)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
5 C u1 p0 c0 {3,S} {4,S} {8,S}
6 C u0 p0 c0 {2,S} {7,D} {9,S}
7 C u1 p0 c0 {4,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-6.58393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67515,0.0530404,-6.24756e-05,3.74708e-08,-8.97632e-12,-709.682,21.0105], Tmin=(100,'K'), Tmax=(1012.19,'K')), NASAPolynomial(coeffs=[10.9787,0.0162757,-7.99455e-06,1.58883e-09,-1.14158e-13,-2593.14,-23.9867], Tmin=(1012.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.58393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s) + radical(CCsJO) + radical(Cds_S)"""),
)

species(
    label = 'F[C]=[C]C1(F)CO1(3099)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
6 C u1 p0 c0 {4,S} {7,D}
7 C u1 p0 c0 {2,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (59.8116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,3150,900,1100,1685,370,167,640,1190,180,180,875.197,875.673,875.984],'cm^-1')),
        HinderedRotor(inertia=(0.191052,'amu*angstrom^2'), symmetry=1, barrier=(4.39266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75678,0.0535265,-7.4526e-05,5.82244e-08,-1.86392e-11,7270.58,21.0056], Tmin=(100,'K'), Tmax=(758.949,'K')), NASAPolynomial(coeffs=[7.8898,0.0212002,-1.0631e-05,2.0942e-09,-1.48289e-13,6339.72,-6.89059], Tmin=(758.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.8116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFH) + ring(Cs-Cs(F)-O2s) + radical(Cds_S) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]C=C1OC1F(3070)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5 C u0 p0 c0 {3,S} {4,S} {6,D}
6 C u0 p0 c0 {5,D} {7,S} {9,S}
7 C u2 p0 c0 {2,S} {6,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (0.098075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33111,0.0409154,-3.58712e-05,1.63609e-08,-3.19287e-12,68.2948,19.0559], Tmin=(100,'K'), Tmax=(1161,'K')), NASAPolynomial(coeffs=[7.82902,0.0219736,-1.13988e-05,2.30855e-09,-1.6698e-13,-1208.33,-8.28904], Tmin=(1161,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.098075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(CsCFHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cs(F)-Cd(Cd)-O2s) + radical(CsCCl_triplet)"""),
)

species(
    label = 'F[C](F)C=C1[CH]O1(2902)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {4,D} {7,S} {8,S}
6 C u1 p0 c0 {3,S} {4,S} {9,S}
7 C u1 p0 c0 {1,S} {2,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-129.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,1000,161,297,490,584,780,1358,532.915,533.57,533.688,533.745,534.09,534.135,534.291],'cm^-1')),
        HinderedRotor(inertia=(0.000592452,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68098,0.0321335,2.58268e-05,-7.17455e-08,3.36574e-11,-15496.3,19.1091], Tmin=(100,'K'), Tmax=(936.696,'K')), NASAPolynomial(coeffs=[20.7678,-0.00135276,2.55177e-06,-4.49444e-10,2.13639e-14,-21178.7,-82.9706], Tmin=(936.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(CsCFFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(methyleneoxirane) + radical(C=CCJO) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
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
    E0 = (-0.546298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (129.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-52.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-4.25417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (300.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (254.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (140.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (158.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (381.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (219.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (109.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (131.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (118.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (103.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C]=CC1(F)[CH]O1(2817)'],
    products = ['C2HF(57)', 'FC1=CO1(1231)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(60.1484,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 60.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]=CC1O[C]1F(2818)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]=CC1(F)[CH]O1(2817)'],
    products = ['FC1=CC2(F)OC12(3086)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;Cdsinglepri_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=CC(F)=C[C]F(3094)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', 'F[C]=CC1=CO1(3095)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(27.5456,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]F(2014)', 'FC1=CO1(1231)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2HF(57)', 'F[C]1[CH]O1(2828)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.001023,'m^3/(mol*s)'), n=2.607, Ea=(13.3553,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', 'FC#CC1(F)[CH]O1(3096)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.29124e+14,'m^3/(mol*s)'), n=-1.78703, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_N-Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]F(2014)', 'F[C]1[CH]O1(2828)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction10',
    reactants = ['HF(38)', 'F[C]=C=C1[CH]O1(3097)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(250.726,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[C]=CC1(F)[CH]O1(2817)'],
    products = ['FC=[C]C1(F)[CH]O1(3098)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_single;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C]=[C]C1(F)CO1(3099)'],
    products = ['F[C]=CC1(F)[CH]O1(2817)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2721.93,'s^-1'), n=2.83, Ea=(138.542,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;Cs_H_out_H/NonDeO] + [R3H_SS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_SS_23cy3;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C]=CC1(F)[CH]O1(2817)'],
    products = ['F[C]C=C1OC1F(3070)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(178.833,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C]=CC1(F)[CH]O1(2817)'],
    products = ['F[C](F)C=C1[CH]O1(2902)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(164.655,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1267',
    isomers = [
        'F[C]=CC1(F)[CH]O1(2817)',
    ],
    reactants = [
        ('C2HF(57)', 'FC1=CO1(1231)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1267',
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

