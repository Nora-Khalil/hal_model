species(
    label = 'O=CC(O)(F)OC(=O)F(9649)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-1133.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,482,664,788,1296,1923,187.417,187.417,187.417],'cm^-1')),
        HinderedRotor(inertia=(0.920149,'amu*angstrom^2'), symmetry=1, barrier=(22.9352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920148,'amu*angstrom^2'), symmetry=1, barrier=(22.9352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920147,'amu*angstrom^2'), symmetry=1, barrier=(22.9352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49309,'amu*angstrom^2'), symmetry=1, barrier=(37.2162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.805435,0.0758444,-0.000100113,6.81286e-08,-1.86435e-11,-136268,26.7759], Tmin=(100,'K'), Tmax=(887.393,'K')), NASAPolynomial(coeffs=[12.2481,0.0242654,-1.29258e-05,2.62771e-09,-1.90217e-13,-138299,-27.0612], Tmin=(887.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1133.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COFOO)"""),
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
    label = 'CO2(14)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'O=C(O)F(635)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(1.70004,'amu*angstrom^2'), symmetry=1, barrier=(39.0874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96655,0.00185137,4.47909e-05,-7.96564e-08,4.23622e-11,-75247.4,7.80392], Tmin=(10,'K'), Tmax=(614.343,'K')), NASAPolynomial(coeffs=[2.83475,0.0173246,-1.27762e-05,4.2862e-09,-5.35304e-13,-75261.2,11.4681], Tmin=(614.343,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-625.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC(O)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)OC(O)F(7845)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-1020.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,482,664,788,1296,1923,222.971,223.169,1069.02],'cm^-1')),
        HinderedRotor(inertia=(0.265891,'amu*angstrom^2'), symmetry=1, barrier=(9.39111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817569,'amu*angstrom^2'), symmetry=1, barrier=(28.8891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608996,'amu*angstrom^2'), symmetry=1, barrier=(21.5108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3648.6,'J/mol'), sigma=(5.56031,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.90 K, Pc=48.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71583,0.0500006,-5.30385e-05,2.7775e-08,-5.77561e-12,-122651,21.0141], Tmin=(100,'K'), Tmax=(1161.05,'K')), NASAPolynomial(coeffs=[12.0233,0.0144896,-7.16047e-06,1.43206e-09,-1.03373e-13,-125045,-30.2527], Tmin=(1161.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1020.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO)"""),
)

species(
    label = 'O=CC(O)(F)OF(8032)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {2,S} {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
7 C u0 p0 c0 {5,D} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-625.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.14501,'amu*angstrom^2'), symmetry=1, barrier=(26.326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14597,'amu*angstrom^2'), symmetry=1, barrier=(26.3481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.146,'amu*angstrom^2'), symmetry=1, barrier=(26.3488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17892,0.0677153,-0.000104804,8.11003e-08,-2.39716e-11,-75164.7,21.2063], Tmin=(100,'K'), Tmax=(675.273,'K')), NASAPolynomial(coeffs=[10.6962,0.0180934,-9.58077e-06,1.9019e-09,-1.34272e-13,-76604,-22.1125], Tmin=(675.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-625.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH)"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787932,2.90876e-06,-1.47517e-09,2.12837e-13,-30281.6,-0.311363], Tmin=(100,'K'), Tmax=(1130.23,'K')), NASAPolynomial(coeffs=[2.84325,0.00275108,-7.81029e-07,1.07243e-10,-5.79388e-15,-29958.6,5.91041], Tmin=(1130.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C=C(F)OC(=O)F(5155)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {8,D}
7 C u0 p0 c0 {2,S} {3,S} {4,D}
8 C u0 p0 c0 {5,D} {6,D}
"""),
    E0 = (-732.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([197,221,431,657,482,664,788,1296,1923,2120,512.5,787.5,180,180,563.499,899.015],'cm^-1')),
        HinderedRotor(inertia=(0.00827737,'amu*angstrom^2'), symmetry=1, barrier=(4.74744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31462,'amu*angstrom^2'), symmetry=1, barrier=(30.2258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.027,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45625,0.0612953,-9.19041e-05,7.04041e-08,-2.15631e-11,-87991.6,21.8663], Tmin=(100,'K'), Tmax=(797.465,'K')), NASAPolynomial(coeffs=[10.1152,0.0178631,-1.021e-05,2.10936e-09,-1.5316e-13,-89372.7,-17.9485], Tmin=(797.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-732.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + group(COFOO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C(F)OC(O)=COF(9785)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {10,S}
9  C u0 p0 c0 {1,S} {3,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-724.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,350,440,435,1725,3010,987.5,1337.5,450,1655,482,664,788,1296,1923,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.853817,'amu*angstrom^2'), symmetry=1, barrier=(19.6309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853688,'amu*angstrom^2'), symmetry=1, barrier=(19.628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853281,'amu*angstrom^2'), symmetry=1, barrier=(19.6186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853468,'amu*angstrom^2'), symmetry=1, barrier=(19.6229,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.255195,0.0934888,-0.000131604,8.64108e-08,-2.16556e-11,-87034.9,27.2055], Tmin=(100,'K'), Tmax=(990.074,'K')), NASAPolynomial(coeffs=[20.4917,0.00966973,-4.61522e-06,9.03711e-10,-6.45975e-14,-91143.1,-72.6794], Tmin=(990.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-724.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(COFOO)"""),
)

species(
    label = 'O=C(F)OOC=C(O)F(9801)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {3,S} {8,D} {10,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u0 p0 c0 {2,S} {4,S} {6,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-861.845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,482,664,788,1296,1923,180],'cm^-1')),
        HinderedRotor(inertia=(1.44138,'amu*angstrom^2'), symmetry=1, barrier=(33.1403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44139,'amu*angstrom^2'), symmetry=1, barrier=(33.1404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44168,'amu*angstrom^2'), symmetry=1, barrier=(33.1471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44108,'amu*angstrom^2'), symmetry=1, barrier=(33.1333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.164618,0.0723302,-6.43916e-05,1.86564e-08,8.61713e-13,-103507,30.3864], Tmin=(100,'K'), Tmax=(1023.86,'K')), NASAPolynomial(coeffs=[22.3788,0.00823347,-3.72793e-06,8.00448e-10,-6.32455e-14,-109246,-83.1159], Tmin=(1023.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-861.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COFOO)"""),
)

species(
    label = 'O=C(F)OC(F)=COO(9802)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-759.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,326,540,652,719,1357,3010,987.5,1337.5,450,1655,482,664,788,1296,1923,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.34878,'amu*angstrom^2'), symmetry=1, barrier=(31.0111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538535,'amu*angstrom^2'), symmetry=1, barrier=(12.382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145371,'amu*angstrom^2'), symmetry=1, barrier=(3.34237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34858,'amu*angstrom^2'), symmetry=1, barrier=(31.0064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502116,0.0845738,-0.000130416,1.06663e-07,-3.52276e-11,-91229.3,28.2989], Tmin=(100,'K'), Tmax=(738.869,'K')), NASAPolynomial(coeffs=[10.9349,0.0280982,-1.57715e-05,3.22911e-09,-2.32883e-13,-92771.1,-18.8768], Tmin=(738.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-759.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(COFOO)"""),
)

species(
    label = 'OC1(F)[CH]OO[C](F)O1(9803)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u1 p0 c0 {4,S} {7,S} {10,S}
9  C u1 p0 c0 {2,S} {3,S} {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-559.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,364,474,579,619,655,1345,2950,1000,482,586,761,1411,538.871,1002.44,1002.44,1002.44,1002.44,1002.44,1002.44,1002.44,1002.44,1002.44,2322.64],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.760386,0.0788591,-8.74647e-05,4.62219e-08,-8.80833e-12,-67039.8,28.4298], Tmin=(100,'K'), Tmax=(1580.72,'K')), NASAPolynomial(coeffs=[19.3388,0.00560599,3.29656e-06,-1.01824e-09,8.03145e-14,-70596.5,-68.8913], Tmin=(1580.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-559.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsFHOO) + ring(124trioxane) + radical(CCsJOOC) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]CC([O])(F)OC(=O)F(9804)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-739.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([373,336,502,477,708,1220,1477,2750,2850,1437.5,1250,1305,750,350,482,664,788,1296,1923,180,180,180,180,2151.61],'cm^-1')),
        HinderedRotor(inertia=(0.230727,'amu*angstrom^2'), symmetry=1, barrier=(5.30486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60295,'amu*angstrom^2'), symmetry=1, barrier=(36.855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0113178,'amu*angstrom^2'), symmetry=1, barrier=(36.9576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677375,0.0831548,-0.000143385,1.35025e-07,-4.97285e-11,-88871.4,29.431], Tmin=(100,'K'), Tmax=(806.13,'K')), NASAPolynomial(coeffs=[6.78668,0.0345062,-1.8747e-05,3.73598e-09,-2.63051e-13,-89260.7,4.9687], Tmin=(806.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-739.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(COFOO) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCOJ)"""),
)

species(
    label = '[O]C(F)([CH]O)OC(=O)F(9805)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u1 p0 c0 {4,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-785.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,3025,407.5,1350,352.5,482,664,788,1296,1923,276.806,276.806,276.806,2337.61],'cm^-1')),
        HinderedRotor(inertia=(0.180803,'amu*angstrom^2'), symmetry=1, barrier=(9.83067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58626,'amu*angstrom^2'), symmetry=1, barrier=(31.8763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180802,'amu*angstrom^2'), symmetry=1, barrier=(9.83067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58626,'amu*angstrom^2'), symmetry=1, barrier=(31.8763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0953383,0.093997,-0.000159685,1.38141e-07,-4.66527e-11,-94309.7,30.212], Tmin=(100,'K'), Tmax=(811.43,'K')), NASAPolynomial(coeffs=[12.0494,0.0249749,-1.34325e-05,2.65002e-09,-1.85019e-13,-95917.4,-22.9139], Tmin=(811.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-785.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(COFOO) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOH)"""),
)

species(
    label = '[O]C(F)OC([O])(F)C=O(9806)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  C u0 p0 c0 {6,D} {7,S} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-668.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,451,553,637,1069,1180,1265,1301,3056,2782.5,750,1395,475,1775,1000,215.256,215.879,216.086,216.244],'cm^-1')),
        HinderedRotor(inertia=(0.504338,'amu*angstrom^2'), symmetry=1, barrier=(17.0247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290034,'amu*angstrom^2'), symmetry=1, barrier=(9.49025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18659,'amu*angstrom^2'), symmetry=1, barrier=(39.723,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.407916,0.0859186,-0.000139011,1.18244e-07,-3.98161e-11,-80224.6,29.8699], Tmin=(100,'K'), Tmax=(791.178,'K')), NASAPolynomial(coeffs=[11.0696,0.0258645,-1.3492e-05,2.6519e-09,-1.85636e-13,-81719.1,-17.8527], Tmin=(791.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-668.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHOO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(O2sj(Cs-F1sO2sH))"""),
)

species(
    label = '[O]C(F)OC(O)(F)[C]=O(9807)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-757.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,451,553,637,1069,1180,1265,1301,3056,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.703895,'amu*angstrom^2'), symmetry=1, barrier=(16.1839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04816,'amu*angstrom^2'), symmetry=1, barrier=(24.0992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76935,'amu*angstrom^2'), symmetry=1, barrier=(40.6808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701508,'amu*angstrom^2'), symmetry=1, barrier=(16.1291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437869,0.0861652,-0.000131524,9.94675e-08,-2.79678e-11,-90938.1,29.6865], Tmin=(100,'K'), Tmax=(651.249,'K')), NASAPolynomial(coeffs=[12.207,0.0242755,-1.29221e-05,2.57163e-09,-1.81859e-13,-92691.5,-23.7383], Tmin=(651.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-757.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHOO) + group(Cds-OdCsH) + radical(O2sj(Cs-F1sO2sH)) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C(F)(C=O)O[C](O)F(9808)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u0 p2 c0 {9,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {10,S}
9  C u1 p0 c0 {2,S} {3,S} {4,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-718.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,482,586,761,1411,246.486,246.488,246.49,246.492],'cm^-1')),
        HinderedRotor(inertia=(0.794202,'amu*angstrom^2'), symmetry=1, barrier=(34.2417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31888,'amu*angstrom^2'), symmetry=1, barrier=(13.7492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318895,'amu*angstrom^2'), symmetry=1, barrier=(13.7492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682197,'amu*angstrom^2'), symmetry=1, barrier=(29.4126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.238376,0.0890708,-0.000140203,1.12438e-07,-3.55924e-11,-86315.6,30.8208], Tmin=(100,'K'), Tmax=(776.779,'K')), NASAPolynomial(coeffs=[13.1351,0.0226586,-1.19564e-05,2.36934e-09,-1.67094e-13,-88319.1,-28.1405], Tmin=(776.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-718.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHOO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = 'O=[C]C(O)(F)O[C](O)F(9809)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u1 p0 c0 {2,S} {3,S} {5,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-807.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,482,586,761,1411,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.134971,0.0912615,-0.000141716,1.09699e-07,-3.33423e-11,-97023.4,31.0949], Tmin=(100,'K'), Tmax=(809.558,'K')), NASAPolynomial(coeffs=[14.3843,0.0208519,-1.12488e-05,2.25434e-09,-1.60294e-13,-99330.4,-34.6387], Tmin=(809.558,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-807.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFHOO) + group(Cds-OdCsH) + radical(Csj(F1s)(O2s-Cs)(O2s-H)) + radical(CsCJ=O)"""),
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
    label = '[O]C=C(O)OC(=O)F(9810)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {5,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-741.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,482,664,788,1296,1923,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02686,'amu*angstrom^2'), symmetry=1, barrier=(23.6095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0274,'amu*angstrom^2'), symmetry=1, barrier=(23.622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02599,'amu*angstrom^2'), symmetry=1, barrier=(23.5895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.130091,0.0797043,-0.000100277,5.67045e-08,-1.19356e-11,-88995.8,25.6942], Tmin=(100,'K'), Tmax=(1260.77,'K')), NASAPolynomial(coeffs=[23.4717,-0.000417044,1.28286e-06,-2.95101e-10,2.06748e-14,-94530.5,-91.9879], Tmin=(1260.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-741.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(COFOO) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]OC(O)(F)C=O(9811)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u0 p0 c0 {4,D} {6,S} {9,S}
8  C u1 p0 c0 {2,S} {5,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-714.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.4525,'amu*angstrom^2'), symmetry=1, barrier=(33.3958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45405,'amu*angstrom^2'), symmetry=1, barrier=(33.4316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4502,'amu*angstrom^2'), symmetry=1, barrier=(33.3429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45188,'amu*angstrom^2'), symmetry=1, barrier=(33.3815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997024,0.0716246,-0.000104253,8.03232e-08,-2.48946e-11,-85874.9,24.0207], Tmin=(100,'K'), Tmax=(787.602,'K')), NASAPolynomial(coeffs=[10.5065,0.0233255,-1.22605e-05,2.45076e-09,-1.74722e-13,-87372.8,-19.5859], Tmin=(787.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-714.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[O]C(=O)F(459)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-381.43,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(62.9882,'amu')),
        NonlinearRotor(inertia=([36.4335,44.7291,81.1626],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([509.533,537.43,765.261,1011.27,1202.3,1550.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C[C](O)F(4284)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)(F)C=O(8194)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-516.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.860907,'amu*angstrom^2'), symmetry=1, barrier=(19.7939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855828,'amu*angstrom^2'), symmetry=1, barrier=(19.6772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73969,0.0545243,-9.06153e-05,7.81645e-08,-2.61901e-11,-62047.4,19.7878], Tmin=(100,'K'), Tmax=(841.881,'K')), NASAPolynomial(coeffs=[8.15902,0.016084,-7.97774e-06,1.52258e-09,-1.04197e-13,-62846.9,-8.40565], Tmin=(841.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-516.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
)

species(
    label = 'OH(7)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92694e-05,-5.32137e-07,1.01946e-09,-3.85933e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07193,0.000604024,-1.3983e-08,-2.13435e-11,2.48057e-15,3579.39,4.57803], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C[C](F)OC(=O)F(4683)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {1,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-786.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,482,664,788,1296,1923,344.877,344.926,344.928],'cm^-1')),
        HinderedRotor(inertia=(0.0309748,'amu*angstrom^2'), symmetry=1, barrier=(52.317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410129,'amu*angstrom^2'), symmetry=1, barrier=(34.6188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.619714,'amu*angstrom^2'), symmetry=1, barrier=(52.3185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12573,0.0457492,-4.28626e-05,2.13068e-08,-4.53568e-12,-94498.3,21.5012], Tmin=(100,'K'), Tmax=(1078.47,'K')), NASAPolynomial(coeffs=[8.01391,0.0239104,-1.24881e-05,2.5306e-09,-1.83208e-13,-95768.3,-7.35075], Tmin=(1078.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-786.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COFOO) + radical(CsCOF1sO2s)"""),
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
    label = '[O]C(F)(C=O)OC(=O)F(9812)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {7,S} {9,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {6,D}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-890.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,482,664,788,1296,1923,280.788,280.788,280.788,280.788],'cm^-1')),
        HinderedRotor(inertia=(0.59185,'amu*angstrom^2'), symmetry=1, barrier=(33.1128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268732,'amu*angstrom^2'), symmetry=1, barrier=(15.035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59185,'amu*angstrom^2'), symmetry=1, barrier=(33.1128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861781,0.0752033,-0.000112453,8.82762e-08,-2.78252e-11,-106957,27.3469], Tmin=(100,'K'), Tmax=(774.984,'K')), NASAPolynomial(coeffs=[10.8316,0.023745,-1.28534e-05,2.5971e-09,-1.86108e-13,-108502,-18.2102], Tmin=(774.984,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COFOO) + radical(C=OCOJ)"""),
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
    label = 'O=C(F)O[C](O)F(1886)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {3,S} {4,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-814.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,586,761,1411,482,664,788,1296,1923,189.176,1054.98,1905.38],'cm^-1')),
        HinderedRotor(inertia=(0.583046,'amu*angstrom^2'), symmetry=1, barrier=(14.7893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08391,'amu*angstrom^2'), symmetry=1, barrier=(27.4997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58374,'amu*angstrom^2'), symmetry=1, barrier=(40.1069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00526,0.0434888,-4.46044e-05,2.23764e-08,-4.46436e-12,-97863,22.1567], Tmin=(100,'K'), Tmax=(1205,'K')), NASAPolynomial(coeffs=[11.2161,0.012913,-6.54267e-06,1.31844e-09,-9.54185e-14,-100083,-23.9976], Tmin=(1205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO) + radical(CsF1sO2sO2s)"""),
)

species(
    label = 'O=[C]C(O)(F)OC(=O)F(9813)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,D}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-979.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,482,664,788,1296,1923,1855,455,950,187.428,187.605,187.89],'cm^-1')),
        HinderedRotor(inertia=(1.05208,'amu*angstrom^2'), symmetry=1, barrier=(24.2662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97719,'amu*angstrom^2'), symmetry=1, barrier=(24.3496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52171,'amu*angstrom^2'), symmetry=1, barrier=(37.3601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.965595,'amu*angstrom^2'), symmetry=1, barrier=(24.1441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763985,0.077334,-0.000113786,8.53479e-08,-2.55118e-11,-117665,27.6005], Tmin=(100,'K'), Tmax=(817.785,'K')), NASAPolynomial(coeffs=[12.0961,0.0219079,-1.21263e-05,2.47713e-09,-1.78873e-13,-119518,-24.7914], Tmin=(817.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-979.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COFOO) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C=C(O)F(1841)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-327.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,197,221,431,657,2120,512.5,787.5,516.311],'cm^-1')),
        HinderedRotor(inertia=(0.000632379,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88591,0.00765106,8.23423e-05,-2.20006e-07,1.65229e-10,-39431.7,9.37981], Tmin=(10,'K'), Tmax=(476.863,'K')), NASAPolynomial(coeffs=[5.30471,0.0166965,-1.1999e-05,3.99878e-09,-4.98819e-13,-39805.2,1.08843], Tmin=(476.863,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-327.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCDC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CC(=O)OC(=O)F(9814)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,D} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u0 p0 c0 {1,S} {2,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-801.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,482,664,788,1296,1923,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33561,0.0642448,-0.000101826,8.52191e-08,-2.86604e-11,-96292.6,25.8333], Tmin=(100,'K'), Tmax=(727.271,'K')), NASAPolynomial(coeffs=[9.30981,0.0203906,-1.13848e-05,2.3215e-09,-1.66757e-13,-97452.6,-10.0988], Tmin=(727.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-801.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + group(COFOO)"""),
)

species(
    label = 'O=C=C(O)OC(=O)F(9815)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,S} {8,D}
7 C u0 p0 c0 {1,S} {2,S} {4,D}
8 C u0 p0 c0 {5,D} {6,D}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-733.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,482,664,788,1296,1923,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3057,'amu*angstrom^2'), symmetry=1, barrier=(30.0205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30077,'amu*angstrom^2'), symmetry=1, barrier=(29.9072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29405,'amu*angstrom^2'), symmetry=1, barrier=(29.7528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71035,0.0858968,-0.000104397,5.28997e-08,-9.4263e-12,-87999.8,29.9055], Tmin=(100,'K'), Tmax=(1697.2,'K')), NASAPolynomial(coeffs=[28.1163,-0.0134995,9.16953e-06,-1.81224e-09,1.20977e-13,-93933,-117.421], Tmin=(1697.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-733.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + group(COFOO) + missing(Cdd-CdO2d)"""),
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
    E0 = (-505.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-297.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-122.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-269.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-51.7115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-212.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-141.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-28.1493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-145.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-229.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-112.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-201.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-162.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-251.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-137.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-111.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-267.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-176.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-226.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-147.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-250.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-233.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-435.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-352.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-407.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-275.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=CC(O)(F)OC(=O)F(9649)'],
    products = ['HF(38)', 'CO2(14)', 'O=CC(=O)F(4234)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.29418e+10,'s^-1'), n=0.745914, Ea=(97.1116,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'O=C(F)OC(O)F(7845)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.81927e-16,'m^3/(mol*s)'), n=6.80628, Ea=(310.857,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CO(13)', 'O=CC(O)(F)OF(8032)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(91.1518,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H2O(3)', 'O=C=C(F)OC(=O)F(5155)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.75748e-06,'m^3/(mol*s)'), n=3.49511, Ea=(183.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OH] for rate rule [Cdd_Cd;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)OC(O)=COF(9785)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(142.352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C(F)OOC=C(O)F(9801)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(117.996,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C(F)OC(F)=COO(9802)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(86.8622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OC1(F)[CH]OO[C](F)O1(9803)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]CC([O])(F)OC(=O)F(9804)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)([CH]O)OC(=O)F(9805)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)OC([O])(F)C=O(9806)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)OC(O)(F)[C]=O(9807)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(C=O)O[C](O)F(9808)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]C(O)(F)O[C](O)F(9809)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C=C(O)OC(=O)F(9810)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.13992e+08,'m^3/(mol*s)'), n=-0.108893, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R_Ext-1C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=[C]OC(O)(F)C=O(9811)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(=O)F(459)', 'O=C[C](O)F(4284)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(51)', '[O]C(O)(F)C=O(8194)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OH(7)', 'O=C[C](F)OC(=O)F(4683)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R
Ea raised from -1.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(F)(C=O)OC(=O)F(9812)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HCO(15)', 'O=C(F)O[C](O)F(1886)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', 'O=[C]C(O)(F)OC(=O)F(9813)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(2.64716,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=CC(O)(F)OC(=O)F(9649)'],
    products = ['O=C(O)F(635)', 'O=CC(=O)F(4234)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.82035e+06,'s^-1'), n=1.67734, Ea=(167.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2285837836968013, var=2.96975919141528, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=CC(O)(F)OC(=O)F(9649)'],
    products = ['O=C(O)F(635)', 'O=C=C(O)F(1841)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.45626e+12,'s^-1'), n=0, Ea=(250.803,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R_1BrClFINOPSSi->O_Ext-2R!H-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R_1BrClFINOPSSi->O_Ext-2R!H-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'O=CC(=O)OC(=O)F(9814)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(144.124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'O=C=C(O)OC(=O)F(9815)'],
    products = ['O=CC(O)(F)OC(=O)F(9649)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(208.516,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #2911',
    isomers = [
        'O=CC(O)(F)OC(=O)F(9649)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'O=CC(=O)F(4234)'),
        ('O=C(O)F(635)', 'O=CC(=O)F(4234)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2911',
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

