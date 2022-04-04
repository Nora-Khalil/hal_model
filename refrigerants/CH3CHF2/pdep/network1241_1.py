species(
    label = 'C#COOCC(F)O[O](3535)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-10.2824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.1358,'amu*angstrom^2'), symmetry=1, barrier=(26.1144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43768,'amu*angstrom^2'), symmetry=1, barrier=(33.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43657,'amu*angstrom^2'), symmetry=1, barrier=(33.0296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43697,'amu*angstrom^2'), symmetry=1, barrier=(33.0387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13693,'amu*angstrom^2'), symmetry=1, barrier=(26.1403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0369651,0.0968741,-0.000138919,9.91188e-08,-2.47275e-11,-1103.11,29.6178], Tmin=(100,'K'), Tmax=(620.307,'K')), NASAPolynomial(coeffs=[11.9929,0.0332083,-1.74431e-05,3.46933e-09,-2.45932e-13,-2844.78,-24.4359], Tmin=(620.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.2824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
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
    label = '[O]OC(F)C=O(1537)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-328.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,232,360,932,1127,1349,1365,3045,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.18317,'amu*angstrom^2'), symmetry=1, barrier=(27.2033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603055,'amu*angstrom^2'), symmetry=1, barrier=(13.8654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3728.91,'J/mol'), sigma=(5.84942,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.45 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15439,0.0440801,-5.80197e-05,4.31461e-08,-1.32601e-11,-39422,18.3048], Tmin=(100,'K'), Tmax=(786.755,'K')), NASAPolynomial(coeffs=[7.18333,0.0185099,-9.26418e-06,1.82889e-09,-1.29953e-13,-40213.3,-4.7504], Tmin=(786.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'CH2(S)(25)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896,2896.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144066,5.45064e-06,-3.57997e-09,7.56176e-13,50400.6,-0.411759], Tmin=(100,'K'), Tmax=(1442.38,'K')), NASAPolynomial(coeffs=[2.62651,0.00394757,-1.49921e-06,2.54533e-10,-1.62951e-14,50691.7,6.78356], Tmin=(1442.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#COOC(F)O[O](3547)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
7  C u0 p0 c0 {4,S} {8,T}
8  C u0 p0 c0 {7,T} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-14.9904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,509,613,660,1171,1360,1414,3084,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.4496,'amu*angstrom^2'), symmetry=1, barrier=(33.3293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45117,'amu*angstrom^2'), symmetry=1, barrier=(33.3654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45088,'amu*angstrom^2'), symmetry=1, barrier=(33.3586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44982,'amu*angstrom^2'), symmetry=1, barrier=(33.3343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624923,0.0827841,-0.000125505,8.89981e-08,-2.06375e-11,-1689.49,25.2305], Tmin=(100,'K'), Tmax=(614.645,'K')), NASAPolynomial(coeffs=[11.9595,0.0227097,-1.2304e-05,2.45018e-09,-1.72837e-13,-3341.41,-26.039], Tmin=(614.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.9904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(CsFHOO) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(11)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.49012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1064.4,1465.7,3224.93],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C#COOC[C]F-2(1858)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
5  C u0 p1 c0 {1,S} {4,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (329.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,617,898,1187,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.41053,'amu*angstrom^2'), symmetry=1, barrier=(32.4309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41111,'amu*angstrom^2'), symmetry=1, barrier=(32.4441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41063,'amu*angstrom^2'), symmetry=1, barrier=(32.433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41022,'amu*angstrom^2'), symmetry=1, barrier=(32.4237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11413,0.0687813,-9.93822e-05,7.7365e-08,-2.43418e-11,39716.4,23.5723], Tmin=(100,'K'), Tmax=(775.133,'K')), NASAPolynomial(coeffs=[9.8368,0.0237658,-1.22646e-05,2.43298e-09,-1.72687e-13,38364.2,-16.2871], Tmin=(775.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(CJ2_singlet-FCs) + group(Ct-CtOs) + group(Ct-CtH)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'C#COOCO[O](3548)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {6,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (220.845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.38002,'amu*angstrom^2'), symmetry=1, barrier=(31.7294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38343,'amu*angstrom^2'), symmetry=1, barrier=(31.8078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37907,'amu*angstrom^2'), symmetry=1, barrier=(31.7074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3806,'amu*angstrom^2'), symmetry=1, barrier=(31.7428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0405612,0.0721344,-8.25593e-05,4.29264e-08,-8.28427e-12,26716.9,26.7733], Tmin=(100,'K'), Tmax=(1422.49,'K')), NASAPolynomial(coeffs=[22.2942,0.0013095,8.22577e-07,-2.27832e-10,1.64216e-14,21220.3,-85.4969], Tmin=(1422.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = 'C#COOC(CF)O[O](3549)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-13.8123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.011691,0.0973067,-0.000141659,1.03851e-07,-2.76378e-11,-1526.62,30.9276], Tmin=(100,'K'), Tmax=(627.351,'K')), NASAPolynomial(coeffs=[12.065,0.0329046,-1.74417e-05,3.4833e-09,-2.47508e-13,-3283.96,-23.5556], Tmin=(627.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.8123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
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
    label = 'C#COOCC([O])F(3550)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-7.21844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.938217,'amu*angstrom^2'), symmetry=1, barrier=(21.5715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74669,'amu*angstrom^2'), symmetry=1, barrier=(40.1598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934341,'amu*angstrom^2'), symmetry=1, barrier=(21.4823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.75186,'amu*angstrom^2'), symmetry=1, barrier=(40.2788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.57085,0.0816094,-0.000114055,8.62732e-08,-2.64505e-11,-750.319,27.0547], Tmin=(100,'K'), Tmax=(793.951,'K')), NASAPolynomial(coeffs=[10.9997,0.0290697,-1.47963e-05,2.92978e-09,-2.08117e-13,-2406.37,-20.8525], Tmin=(793.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.21844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(CsCFHO) + group(Ct-CtOs) + group(Ct-CtH) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = 'HCCO(21)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (166.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87608,0.0221205,-3.58869e-05,3.05403e-08,-1.01281e-11,20163.4,13.6968], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91479,0.00371409,-1.30137e-06,2.06473e-10,-1.21477e-14,19359.6,-5.50567], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(166.705,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'FC1COOO1(482)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-275.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27847,0.0247624,2.97527e-05,-6.8331e-08,3.25679e-11,-33006.5,13.5316], Tmin=(100,'K'), Tmax=(880.085,'K')), NASAPolynomial(coeffs=[14.677,0.00466565,2.21309e-06,-6.61897e-10,4.91985e-14,-36592.9,-52.6774], Tmin=(880.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(Cs-CsOsHH) + ring(123trioxolane)"""),
)

species(
    label = 'C#COOC=CF(1855)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,D} {8,S}
5  C u0 p0 c0 {1,S} {4,D} {9,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {6,T} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (129.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.52592,'amu*angstrom^2'), symmetry=1, barrier=(35.084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52882,'amu*angstrom^2'), symmetry=1, barrier=(35.1505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52899,'amu*angstrom^2'), symmetry=1, barrier=(35.1546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5194,0.0590392,-6.93071e-05,4.32907e-08,-1.11155e-11,15680.6,22.0182], Tmin=(100,'K'), Tmax=(933.568,'K')), NASAPolynomial(coeffs=[9.82756,0.0234405,-1.21074e-05,2.44264e-09,-1.76489e-13,14129.4,-17.4925], Tmin=(933.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'FC1COO[C]=COO1(3551)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {3,S} {9,S}
6  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {13,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-44.3016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40672,0.0404211,2.36695e-05,-5.69033e-08,2.29088e-11,-5220.49,25.5617], Tmin=(100,'K'), Tmax=(1037.79,'K')), NASAPolynomial(coeffs=[15.4596,0.0256047,-1.1788e-05,2.40888e-09,-1.80328e-13,-10256.2,-52.9652], Tmin=(1037.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.3016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1OOCC(F)OO1(3552)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {3,S} {8,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u1 p0 c0 {8,D} {13,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-0.924725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666365,0.0326817,0.000106987,-1.96653e-07,8.72957e-11,46.6821,27.8761], Tmin=(100,'K'), Tmax=(917.631,'K')), NASAPolynomial(coeffs=[35.9567,-0.0153574,1.25795e-05,-2.42581e-09,1.51283e-13,-10884.2,-163.616], Tmin=(917.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.924725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(Cds_P)"""),
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
    label = 'C#COOC[CH]O[O](3553)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {5,S} {11,S}
7  C u0 p0 c0 {2,S} {8,T}
8  C u0 p0 c0 {7,T} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (404.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(2.19631,'amu*angstrom^2'), symmetry=1, barrier=(50.4975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19639,'amu*angstrom^2'), symmetry=1, barrier=(50.4993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.492202,'amu*angstrom^2'), symmetry=1, barrier=(11.3167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992835,'amu*angstrom^2'), symmetry=1, barrier=(22.8272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19495,'amu*angstrom^2'), symmetry=1, barrier=(50.4663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0807459,0.0966785,-0.000166565,1.49932e-07,-5.21117e-11,48722.1,29.3951], Tmin=(100,'K'), Tmax=(841.598,'K')), NASAPolynomial(coeffs=[9.44156,0.0323789,-1.66569e-05,3.21571e-09,-2.21223e-13,47848.1,-9.98331], Tmin=(841.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[O]CC(F)O[O](2005)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-192.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.151868,'amu*angstrom^2'), symmetry=1, barrier=(3.49174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151696,'amu*angstrom^2'), symmetry=1, barrier=(3.48778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88506,0.0547729,-9.52024e-05,9.34081e-08,-3.51294e-11,-23131.5,20.6619], Tmin=(100,'K'), Tmax=(835.304,'K')), NASAPolynomial(coeffs=[3.6698,0.0280021,-1.44023e-05,2.80119e-09,-1.94137e-13,-22793.9,16.1785], Tmin=(835.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'C#CO[O](374)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {3,T} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (341.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,207.447,207.514],'cm^-1')),
        HinderedRotor(inertia=(1.42755,'amu*angstrom^2'), symmetry=1, barrier=(43.6061,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8409,0.0270256,-4.16772e-05,3.32116e-08,-1.02958e-11,41097.9,11.7897], Tmin=(100,'K'), Tmax=(872.92,'K')), NASAPolynomial(coeffs=[6.69659,0.00695035,-3.04403e-06,5.47509e-10,-3.61747e-14,40516.5,-5.76235], Tmin=(872.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-OsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(F)O[O](222)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-48.0821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.485905,'amu*angstrom^2'), symmetry=1, barrier=(11.1719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486073,'amu*angstrom^2'), symmetry=1, barrier=(11.1758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3366.65,'J/mol'), sigma=(5.62898,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.86 K, Pc=42.83 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33352,0.0400211,-5.12825e-05,3.75379e-08,-1.13901e-11,-5725.96,17.2582], Tmin=(100,'K'), Tmax=(795.118,'K')), NASAPolynomial(coeffs=[6.83717,0.0173642,-8.53954e-06,1.69956e-09,-1.21736e-13,-6442.13,-3.43672], Tmin=(795.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.0012157,5.31615e-06,-4.8944e-09,1.45844e-12,-1038.59,4.68369], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15383,0.00167803,-7.69968e-07,1.51274e-10,-1.08781e-14,-1040.82,6.16752], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#COOC[CH]F(1689)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
5  C u1 p0 c0 {1,S} {4,S} {10,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {6,T} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (156.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.34765,'amu*angstrom^2'), symmetry=1, barrier=(30.9852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34731,'amu*angstrom^2'), symmetry=1, barrier=(30.9772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34771,'amu*angstrom^2'), symmetry=1, barrier=(30.9864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34691,'amu*angstrom^2'), symmetry=1, barrier=(30.968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.072,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.62,'J/mol'), sigma=(5.46284,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.12 K, Pc=44.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671594,0.0795213,-0.000123828,1.04364e-07,-3.53054e-11,18973.6,24.8109], Tmin=(100,'K'), Tmax=(766.274,'K')), NASAPolynomial(coeffs=[10.0768,0.0269011,-1.39233e-05,2.74368e-09,-1.92999e-13,17635.7,-17.3847], Tmin=(766.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(CsCsFHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = 'C2H(22)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]OCC(F)O[O](3554)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u1 p2 c0 {2,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-195.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,180],'cm^-1')),
        HinderedRotor(inertia=(0.428595,'amu*angstrom^2'), symmetry=1, barrier=(9.85425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.814711,'amu*angstrom^2'), symmetry=1, barrier=(18.7318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427141,'amu*angstrom^2'), symmetry=1, barrier=(9.82082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11096,0.0708368,-0.000118587,1.05575e-07,-3.6481e-11,-23366.7,24.1019], Tmin=(100,'K'), Tmax=(841.818,'K')), NASAPolynomial(coeffs=[7.96252,0.0247536,-1.23699e-05,2.36918e-09,-1.62539e-13,-24040.9,-4.92625], Tmin=(841.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]O[CH]F(283)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-10.1775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,580,1155,1237,1373,3147],'cm^-1')),
        HinderedRotor(inertia=(0.627597,'amu*angstrom^2'), symmetry=1, barrier=(14.4297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87623,0.0275715,-4.75678e-05,4.25822e-08,-1.46837e-11,-1186.36,12.8042], Tmin=(100,'K'), Tmax=(836.137,'K')), NASAPolynomial(coeffs=[5.84968,0.00835189,-4.1278e-06,8.02249e-10,-5.55883e-14,-1508.99,0.0352551], Tmin=(836.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.1775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'C#COO[CH2](901)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u1 p0 c0 {1,S} {6,S} {7,S}
4 C u0 p0 c0 {2,S} {5,T}
5 C u0 p0 c0 {4,T} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (382.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.81214,'amu*angstrom^2'), symmetry=1, barrier=(41.6646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81525,'amu*angstrom^2'), symmetry=1, barrier=(41.7361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22571,0.041128,-4.32669e-05,2.39862e-08,-5.46916e-12,46113.1,15.1973], Tmin=(100,'K'), Tmax=(1043.81,'K')), NASAPolynomial(coeffs=[8.57818,0.0167845,-8.28393e-06,1.64292e-09,-1.17754e-13,44786.9,-15.722], Tmin=(1043.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-OsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CsJOOC)"""),
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
    label = 'C#COO[CH]C(F)O[O](3555)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
7  C u1 p0 c0 {2,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (176.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,487,638,688,1119,1325,1387,3149,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.50254,'amu*angstrom^2'), symmetry=1, barrier=(34.5463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50359,'amu*angstrom^2'), symmetry=1, barrier=(34.5705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50171,'amu*angstrom^2'), symmetry=1, barrier=(34.5272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5028,'amu*angstrom^2'), symmetry=1, barrier=(34.5524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50213,'amu*angstrom^2'), symmetry=1, barrier=(34.5369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.265525,0.10315,-0.000173111,1.51409e-07,-5.20209e-11,21331.6,30.9539], Tmin=(100,'K'), Tmax=(804.417,'K')), NASAPolynomial(coeffs=[11.764,0.03133,-1.68053e-05,3.32061e-09,-2.32489e-13,19784.6,-22.0494], Tmin=(804.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(CCsJOOC)"""),
)

species(
    label = 'C#COOC[C](F)O[O](3556)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u1 p2 c0 {4,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {4,S} {6,S}
8  C u0 p0 c0 {3,S} {9,T}
9  C u0 p0 c0 {8,T} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (184.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,395,473,707,1436,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.734483,'amu*angstrom^2'), symmetry=1, barrier=(16.8872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61601,'amu*angstrom^2'), symmetry=1, barrier=(37.1551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61444,'amu*angstrom^2'), symmetry=1, barrier=(37.1192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.733425,'amu*angstrom^2'), symmetry=1, barrier=(16.8629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61599,'amu*angstrom^2'), symmetry=1, barrier=(37.1548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.406708,0.109828,-0.000199331,1.82215e-07,-6.36165e-11,22309,31.9673], Tmin=(100,'K'), Tmax=(845.007,'K')), NASAPolynomial(coeffs=[10.5101,0.0328304,-1.77019e-05,3.45781e-09,-2.38517e-13,21368.1,-13.5124], Tmin=(845.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[C]#COOCC(F)O[O](3557)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u1 p2 c0 {3,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u1 p0 c0 {8,T}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (326.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.26168,'amu*angstrom^2'), symmetry=1, barrier=(52.0005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.699044,'amu*angstrom^2'), symmetry=1, barrier=(16.0724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26349,'amu*angstrom^2'), symmetry=1, barrier=(52.0421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745394,'amu*angstrom^2'), symmetry=1, barrier=(17.1381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26275,'amu*angstrom^2'), symmetry=1, barrier=(52.025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.298876,0.105411,-0.000183412,1.63169e-07,-5.59166e-11,39456.7,31.5209], Tmin=(100,'K'), Tmax=(844.032,'K')), NASAPolynomial(coeffs=[11.1129,0.0313582,-1.63139e-05,3.1515e-09,-2.16487e-13,38241.7,-17.3854], Tmin=(844.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(Acetyl)"""),
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
    label = 'C#COOC=CO[O](3558)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,D} {9,S}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  C u0 p0 c0 {2,S} {8,T}
8  C u0 p0 c0 {7,T} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (362.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.13461,'amu*angstrom^2'), symmetry=1, barrier=(26.0869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13491,'amu*angstrom^2'), symmetry=1, barrier=(26.0939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13548,'amu*angstrom^2'), symmetry=1, barrier=(26.1069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13572,'amu*angstrom^2'), symmetry=1, barrier=(26.1125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478773,0.0776617,-0.000100674,6.38938e-08,-1.57383e-11,43709.6,28.9422], Tmin=(100,'K'), Tmax=(1000.37,'K')), NASAPolynomial(coeffs=[16.1512,0.0149956,-6.70996e-06,1.27478e-09,-8.94373e-14,40574,-46.6739], Tmin=(1000.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
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
    label = 'C#COOCC(=O)F(3559)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,D} {5,S}
7  C u0 p0 c0 {3,S} {8,T}
8  C u0 p0 c0 {7,T} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-168.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (118.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02555,0.0702279,-8.42228e-05,5.2448e-08,-1.32638e-11,-20203.3,24.7832], Tmin=(100,'K'), Tmax=(952.413,'K')), NASAPolynomial(coeffs=[11.7802,0.0250587,-1.3082e-05,2.64996e-09,-1.91908e-13,-22251.8,-26.5772], Tmin=(952.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'C#COO[CH]C(F)OO(3560)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {2,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u1 p0 c0 {3,S} {6,S} {11,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (24.1536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.360856,0.104548,-0.000159723,1.2934e-07,-4.21267e-11,3053.99,30.5506], Tmin=(100,'K'), Tmax=(750.532,'K')), NASAPolynomial(coeffs=[12.9144,0.0337852,-1.82762e-05,3.67854e-09,-2.62583e-13,1061.6,-29.6833], Tmin=(750.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.1536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(CsCFHO) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOOC)"""),
)

species(
    label = '[C]#COOCC(F)OO(3561)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {3,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u1 p0 c0 {8,T}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (174.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,261,493,600,1152,1365,1422,3097,2175,525,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.522658,0.10852,-0.000177033,1.51866e-07,-5.14756e-11,21184.6,31.5663], Tmin=(100,'K'), Tmax=(798.242,'K')), NASAPolynomial(coeffs=[12.5714,0.0332493,-1.74422e-05,3.42536e-09,-2.39409e-13,19401.8,-26.727], Tmin=(798.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFHO) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = 'C#COOC[CH]OOF(3562)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {1,S} {3,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {3,S} {6,S} {12,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {8,T} {13,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (346.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.6442,'amu*angstrom^2'), symmetry=1, barrier=(37.8034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64415,'amu*angstrom^2'), symmetry=1, barrier=(37.8023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64437,'amu*angstrom^2'), symmetry=1, barrier=(37.8072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64573,'amu*angstrom^2'), symmetry=1, barrier=(37.8386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.497298,0.108994,-0.000181645,1.59735e-07,-5.53301e-11,41786.9,32.166], Tmin=(100,'K'), Tmax=(802.237,'K')), NASAPolynomial(coeffs=[11.6274,0.0352203,-1.88005e-05,3.7141e-09,-2.60236e-13,40270.1,-20.9857], Tmin=(802.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCt) + group(O2sFO) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOO)"""),
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
    E0 = (-103.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (360.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (263.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (372.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (175.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (123.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-68.2353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (22.7419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-57.6252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (54.1467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (364.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-104.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (180.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (35.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (249.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (259.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (275.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (283.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (425.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (188.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (34.5233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-5.89369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (203.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (321.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['CH2CO(28)', '[O]OC(F)C=O(1537)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(702.966,'s^-1'), n=2.72887, Ea=(19.5938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(25)', 'C#COOC(F)O[O](3547)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.26474e-05,'m^3/(mol*s)'), n=2.95311, Ea=(69.4891,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_N-2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_N-2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['HO2(11)', 'C#COOC[C]F-2(1858)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(44.7096,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(40)', 'C#COOCO[O](3548)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.12553e-07,'m^3/(mol*s)'), n=3.34134, Ea=(125.462,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['C#COOC(CF)O[O](3549)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', 'C#COOCC([O])F(3550)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['HCCO(21)', 'FC1COOO1(482)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['HO2(11)', 'C#COOC=CF(1855)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.00406e+10,'s^-1'), n=0.563333, Ea=(145.788,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['FC1COO[C]=COO1(3551)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.28217e+09,'s^-1'), n=0.356567, Ea=(65.4205,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;triplebond_intra_H;radadd_intra] for rate rule [R8_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['[CH]=C1OOCC(F)OO1(3552)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.05e+09,'s^-1'), n=0.155, Ea=(177.192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;triplebond_intra_H;radadd_intra] for rate rule [R8;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'C#COOC[CH]O[O](3553)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['HCCO(21)', '[O]CC(F)O[O](2005)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(34.1892,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CO[O](374)', '[CH2]C(F)O[O](222)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)', 'C#COOC[CH]F(1689)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C2H(22)', '[O]OCC(F)O[O](3554)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(99063.5,'m^3/(mol*s)'), n=0.290298, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.10559855350330952, var=1.0655351314124175, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]O[CH]F(283)', 'C#COO[CH2](901)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'C#COO[CH]C(F)O[O](3555)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'C#COOC[C](F)O[O](3556)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[C]#COOCC(F)O[O](3557)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'C#COOC=CO[O](3558)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(219.876,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['OH(7)', 'C#COOCC(=O)F(3559)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#COOCC(F)O[O](3535)'],
    products = ['C#COO[CH]C(F)OO(3560)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.94e+12,'s^-1'), n=-0.24, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 291 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_OOH/H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]#COOCC(F)OO(3561)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.27e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H;Y_rad_out;XH_out] for rate rule [R8H;Ct_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#COOC[CH]OOF(3562)'],
    products = ['C#COOCC(F)O[O](3535)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(88.2866,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1241',
    isomers = [
        'C#COOCC(F)O[O](3535)',
    ],
    reactants = [
        ('CH2CO(28)', '[O]OC(F)C=O(1537)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1241',
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

