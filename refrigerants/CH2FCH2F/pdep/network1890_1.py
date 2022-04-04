species(
    label = '[O]C(C[C]=O)C(F)F(4691)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-419.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,235,523,627,1123,1142,1372,1406,3097,1855,455,950,338.59,339.19,339.863],'cm^-1')),
        HinderedRotor(inertia=(0.101218,'amu*angstrom^2'), symmetry=1, barrier=(8.26416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10129,'amu*angstrom^2'), symmetry=1, barrier=(8.26233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101155,'amu*angstrom^2'), symmetry=1, barrier=(8.26511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661388,0.07907,-0.000118029,9.66391e-08,-3.18125e-11,-50350.6,28.3317], Tmin=(100,'K'), Tmax=(778.99,'K')), NASAPolynomial(coeffs=[10.1772,0.02747,-1.33982e-05,2.58399e-09,-1.79798e-13,-51750.1,-14.667], Tmin=(778.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-419.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCCJ=O)"""),
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
    label = 'O=CC(F)F(990)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-557.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([195,270,1147,1130,1359,1388,1409,3075,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.595284,'amu*angstrom^2'), symmetry=1, barrier=(13.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93683,0.00748989,0.000222302,-1.39306e-06,2.75418e-09,-67109.3,9.32225], Tmin=(10,'K'), Tmax=(169.939,'K')), NASAPolynomial(coeffs=[3.97443,0.0203448,-1.24424e-05,3.60772e-09,-4.02215e-13,-67130.4,8.62375], Tmin=(169.939,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-557.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODCC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF2(43)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = '[O]CC[C]=O(610)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5 C u1 p0 c0 {2,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (22.5151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,180,3016.62],'cm^-1')),
        HinderedRotor(inertia=(0.221227,'amu*angstrom^2'), symmetry=1, barrier=(5.08643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220149,'amu*angstrom^2'), symmetry=1, barrier=(5.06165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3708.82,'J/mol'), sigma=(6.05458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=579.31 K, Pc=37.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42016,0.0433568,-7.64097e-05,7.96031e-08,-3.08672e-11,2756.53,19.4434], Tmin=(100,'K'), Tmax=(859.035,'K')), NASAPolynomial(coeffs=[1.06986,0.0286198,-1.39648e-05,2.65109e-09,-1.80814e-13,3764.27,30.2679], Tmin=(859.035,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.5151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCCJ=O)"""),
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
    label = '[O]C(F)C[C]=O(1422)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-202.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.232125,'amu*angstrom^2'), symmetry=1, barrier=(5.33701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231102,'amu*angstrom^2'), symmetry=1, barrier=(5.31348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3725.89,'J/mol'), sigma=(5.9714,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.98 K, Pc=39.71 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92641,0.0498515,-7.49997e-05,6.41478e-08,-2.18661e-11,-24263.6,21.4294], Tmin=(100,'K'), Tmax=(827.75,'K')), NASAPolynomial(coeffs=[6.65959,0.0203662,-9.585e-06,1.81176e-09,-1.24215e-13,-24820.6,0.857734], Tmin=(827.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCFHO) + group(Cds-OdCsH) + radical(O2sj(Cs-CsF1sH)) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C([O])C(F)F(4687)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {11,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-446.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499588,0.0708037,-7.59012e-05,3.97787e-08,-8.08256e-12,-53517.6,29.0023], Tmin=(100,'K'), Tmax=(1209.16,'K')), NASAPolynomial(coeffs=[17.6229,0.0141582,-5.63034e-06,1.03491e-09,-7.20382e-14,-57658.6,-56.8597], Tmin=(1209.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-446.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,29230.2,5.12616], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=[C]C[CH]C(F)F(833)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u1 p0 c0 {4,S} {5,S} {11,S}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-277.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.203975,'amu*angstrom^2'), symmetry=1, barrier=(4.68978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204059,'amu*angstrom^2'), symmetry=1, barrier=(4.69172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204206,'amu*angstrom^2'), symmetry=1, barrier=(4.69509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.071,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49714,0.0592632,-7.66316e-05,5.68604e-08,-1.74325e-11,-33299.1,25.6696], Tmin=(100,'K'), Tmax=(789.346,'K')), NASAPolynomial(coeffs=[8.20031,0.0252943,-1.20791e-05,2.33942e-09,-1.64451e-13,-34357.3,-5.08349], Tmin=(789.346,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[C]=O(129)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08918,0.00200392,-1.61651e-05,2.55044e-08,-1.16417e-11,52802.7,4.52499], Tmin=(100,'K'), Tmax=(856.118,'K')), NASAPolynomial(coeffs=[0.961586,0.00569052,-3.48048e-06,7.19212e-10,-5.0805e-14,53738.7,21.4665], Tmin=(856.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])C(F)F(4737)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6  C u1 p0 c0 {4,S} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-261.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.644715,'amu*angstrom^2'), symmetry=1, barrier=(14.8233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644714,'amu*angstrom^2'), symmetry=1, barrier=(14.8232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27707,0.0629161,-8.38103e-05,5.75043e-08,-1.56055e-11,-31389.3,18.6937], Tmin=(100,'K'), Tmax=(902.392,'K')), NASAPolynomial(coeffs=[11.5504,0.0173775,-8.11289e-06,1.58013e-09,-1.11974e-13,-33243.4,-29.8139], Tmin=(902.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'O=C1CC(C(F)F)O1(4695)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
8  C u0 p0 c0 {3,S} {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-743.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74985,0.0369581,1.57182e-05,-5.0174e-08,2.31719e-11,-89297.3,22.5129], Tmin=(100,'K'), Tmax=(936.885,'K')), NASAPolynomial(coeffs=[13.4401,0.0186335,-5.51509e-06,9.20926e-10,-6.4914e-14,-92874.1,-40.5221], Tmin=(936.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-743.229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone)"""),
)

species(
    label = 'O=CCC(=O)C(F)F(4835)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u0 p0 c0 {4,D} {5,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-698.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31395,0.0644247,-8.14619e-05,6.15534e-08,-1.96802e-11,-83879.2,25.9664], Tmin=(100,'K'), Tmax=(751.304,'K')), NASAPolynomial(coeffs=[7.51555,0.0314077,-1.5544e-05,3.06287e-09,-2.17592e-13,-84811.1,-2.17951], Tmin=(751.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + group(Cds-OdCsH)"""),
)

species(
    label = 'O=C=CC(O)C(F)F(4836)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {12,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-692.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.334812,0.0876346,-0.000141629,1.19678e-07,-3.93833e-11,-83208,26.0173], Tmin=(100,'K'), Tmax=(848.563,'K')), NASAPolynomial(coeffs=[11.0285,0.0258607,-1.23411e-05,2.31907e-09,-1.57405e-13,-84613.6,-21.4064], Tmin=(848.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-692.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(CsCsFFH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.56331e-09,3.13594e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88611e-07,1.21037e-10,-7.84012e-15,-14180.9,6.71043], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][C]=O(338)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.3074e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'CHF2(82)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]CC=O(1325)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,D} {3,S} {8,S}
5 C u1 p0 c0 {2,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-112.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.440593,'amu*angstrom^2'), symmetry=1, barrier=(10.1301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440129,'amu*angstrom^2'), symmetry=1, barrier=(10.1194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34836,0.0397447,-5.83745e-05,4.996e-08,-1.7121e-11,-13513.9,16.871], Tmin=(100,'K'), Tmax=(829.98,'K')), NASAPolynomial(coeffs=[5.8159,0.0174819,-8.10675e-06,1.52463e-09,-1.0431e-13,-13898.3,1.9402], Tmin=(829.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCCJ=O)"""),
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
    label = 'O=[C]CC(=O)C(F)F(4837)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u1 p0 c0 {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-538.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,195,270,1147,1130,1359,1388,1409,3075,375,552.5,462.5,1710,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.353578,'amu*angstrom^2'), symmetry=1, barrier=(8.12946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353588,'amu*angstrom^2'), symmetry=1, barrier=(8.12968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353587,'amu*angstrom^2'), symmetry=1, barrier=(8.12967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753084,0.0799955,-0.000136453,1.22405e-07,-4.24011e-11,-64615.6,26.5676], Tmin=(100,'K'), Tmax=(845.113,'K')), NASAPolynomial(coeffs=[8.4375,0.0272331,-1.37109e-05,2.63037e-09,-1.8039e-13,-65329.1,-5.74898], Tmin=(845.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCCJ=O)"""),
)

species(
    label = '[O][CH]C(F)F(1636)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-255.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.637092,'amu*angstrom^2'), symmetry=1, barrier=(14.648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0921,0.0467254,-7.94794e-05,7.13576e-08,-2.49602e-11,-30696.7,15.5487], Tmin=(100,'K'), Tmax=(816.59,'K')), NASAPolynomial(coeffs=[6.94864,0.0153515,-7.91621e-06,1.55898e-09,-1.09043e-13,-31237,-5.34886], Tmin=(816.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C(C=C=O)C(F)F(4838)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-462.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.749033,'amu*angstrom^2'), symmetry=1, barrier=(17.2217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745905,'amu*angstrom^2'), symmetry=1, barrier=(17.1498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499274,0.083625,-0.000135695,1.14669e-07,-3.78308e-11,-55507.5,25.235], Tmin=(100,'K'), Tmax=(837.777,'K')), NASAPolynomial(coeffs=[10.9261,0.0241328,-1.17934e-05,2.24073e-09,-1.53186e-13,-56913.9,-21.1891], Tmin=(837.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-462.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(CsCsFFH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CC(C)OJ)"""),
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
    label = 'O=[C]CC(=O)[CH]F(4839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u1 p0 c0 {1,S} {5,S} {10,S}
7  C u1 p0 c0 {3,D} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-177.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,235,1215,1347,1486,3221,1855,455,950,260.017,261.253],'cm^-1')),
        HinderedRotor(inertia=(0.101347,'amu*angstrom^2'), symmetry=1, barrier=(5.32668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109162,'amu*angstrom^2'), symmetry=1, barrier=(5.29262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344821,'amu*angstrom^2'), symmetry=1, barrier=(17.2667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28868,0.0654276,-0.000104453,9.01444e-08,-3.05273e-11,-21219.8,23.0276], Tmin=(100,'K'), Tmax=(836.553,'K')), NASAPolynomial(coeffs=[8.19003,0.0228909,-1.10801e-05,2.10454e-09,-1.44121e-13,-22040.7,-7.04109], Tmin=(836.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(Csj(CO-CsO2d)(F1s)(H)) + radical(C=OCCJ=O)"""),
)

species(
    label = 'O=[C]C[C](O)C(F)F(4840)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {7,S} {12,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,D} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-473.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.274712,0.0926627,-0.000164177,1.48183e-07,-5.08568e-11,-56804.1,29.3327], Tmin=(100,'K'), Tmax=(866.516,'K')), NASAPolynomial(coeffs=[9.22374,0.0287098,-1.42743e-05,2.69736e-09,-1.82246e-13,-57505,-7.65403], Tmin=(866.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-473.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=CC([O])C(F)F(4841)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {5,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {4,S} {7,D} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-436.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519912,0.0647382,-4.67515e-05,3.0752e-09,6.17105e-12,-52379.7,28.6387], Tmin=(100,'K'), Tmax=(982.882,'K')), NASAPolynomial(coeffs=[20.1864,0.0106236,-3.72512e-06,7.23363e-10,-5.55748e-14,-57497.7,-72.2704], Tmin=(982.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-436.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(CsCsFFH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=[C][CH]C(O)C(F)F(4842)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {12,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u1 p0 c0 {5,S} {8,S} {11,S}
8  C u1 p0 c0 {4,D} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-450.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849234,0.0755867,-0.00010491,7.38195e-08,-1.87071e-11,-54021.1,30.3701], Tmin=(100,'K'), Tmax=(656.704,'K')), NASAPolynomial(coeffs=[10.7549,0.0250823,-1.20069e-05,2.30368e-09,-1.60169e-13,-55534.1,-14.8674], Tmin=(656.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-450.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC(O)[C](F)F(4843)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {12,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
"""),
    E0 = (-448.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1855,455,950,225.867,225.867],'cm^-1')),
        HinderedRotor(inertia=(0.176633,'amu*angstrom^2'), symmetry=1, barrier=(6.39448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176633,'amu*angstrom^2'), symmetry=1, barrier=(6.39448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176633,'amu*angstrom^2'), symmetry=1, barrier=(6.39448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176633,'amu*angstrom^2'), symmetry=1, barrier=(6.39448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433861,0.0871587,-0.000147269,1.30018e-07,-4.43508e-11,-53772.2,30.385], Tmin=(100,'K'), Tmax=(848.849,'K')), NASAPolynomial(coeffs=[9.54374,0.0279464,-1.38591e-05,2.64103e-09,-1.80361e-13,-54732.1,-8.61638], Tmin=(848.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-448.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s)) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C](CC=O)C(F)F(4844)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u0 p0 c0 {4,D} {5,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-402.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,522,611,926,1093,1137,1374,1416,3112,360,370,350,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.254061,'amu*angstrom^2'), symmetry=1, barrier=(5.84137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254912,'amu*angstrom^2'), symmetry=1, barrier=(5.86093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439624,'amu*angstrom^2'), symmetry=1, barrier=(10.1078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.52721,0.0855308,-0.000142811,1.28779e-07,-4.5262e-11,-48344.6,27.4148], Tmin=(100,'K'), Tmax=(831.325,'K')), NASAPolynomial(coeffs=[8.16419,0.0320919,-1.62685e-05,3.14638e-09,-2.17464e-13,-49037.5,-4.5491], Tmin=(831.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[O]C(CC=O)[C](F)F(4845)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  C u0 p0 c0 {4,D} {6,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-377.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,266.758,267.188,267.297],'cm^-1')),
        HinderedRotor(inertia=(0.553617,'amu*angstrom^2'), symmetry=1, barrier=(28.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186515,'amu*angstrom^2'), symmetry=1, barrier=(9.44239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00182908,'amu*angstrom^2'), symmetry=1, barrier=(9.44078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702259,0.0798344,-0.000125225,1.09781e-07,-3.84538e-11,-45313.3,28.4105], Tmin=(100,'K'), Tmax=(797.805,'K')), NASAPolynomial(coeffs=[8.39717,0.0314828,-1.5945e-05,3.11221e-09,-2.17447e-13,-46230.2,-5.02592], Tmin=(797.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'O=[C]CC([CH]F)OF(4846)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {5,S} {12,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-99.3849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.214694,'amu*angstrom^2'), symmetry=1, barrier=(4.93623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210971,'amu*angstrom^2'), symmetry=1, barrier=(4.85065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209864,'amu*angstrom^2'), symmetry=1, barrier=(4.82519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781752,'amu*angstrom^2'), symmetry=1, barrier=(17.974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24143,0.0942133,-0.000168305,1.54704e-07,-5.42228e-11,-11829,29.7373], Tmin=(100,'K'), Tmax=(853.534,'K')), NASAPolynomial(coeffs=[8.50041,0.0316011,-1.62555e-05,3.12663e-09,-2.13951e-13,-12368,-3.69811], Tmin=(853.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.3849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(Cds-OdCsH) + radical(Csj(Cs-O2sCsH)(F1s)(H)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C([CH]F)CC(=O)F(4847)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {2,S} {5,S} {12,S}
8  C u0 p0 c0 {1,S} {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-417.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,486,617,768,1157,1926,223.517,224.042,224.244,2435.94],'cm^-1')),
        HinderedRotor(inertia=(0.179264,'amu*angstrom^2'), symmetry=1, barrier=(6.47581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181839,'amu*angstrom^2'), symmetry=1, barrier=(6.48067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978714,'amu*angstrom^2'), symmetry=1, barrier=(34.7732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576966,0.0852141,-0.00014451,1.33052e-07,-4.75959e-11,-50098.2,28.275], Tmin=(100,'K'), Tmax=(830.116,'K')), NASAPolynomial(coeffs=[7.28126,0.0338933,-1.74142e-05,3.38597e-09,-2.34685e-13,-50556.1,1.12516], Tmin=(830.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(COCsFO) + radical(CC(C)OJ) + radical(Csj(Cs-O2sCsH)(F1s)(H))"""),
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
    E0 = (-192.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (172.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (313.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (52.6047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (192.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (403.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-184.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-129.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-129.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-126.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-124.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-105.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-82.7987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-79.7559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-22.7919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (131.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (14.2497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-27.3258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-33.7259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-117.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-192.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-71.2725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-77.5264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (225.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (50.0951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['CH2CO(28)', 'O=CC(F)F(990)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', '[O]CC[C]=O(610)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33604e-10,'m^3/(mol*s)'), n=4.72997, Ea=(126.614,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_3Br1sCCl1sF1sHI1s->F1s_2Br1sCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', '[O]C(F)C[C]=O(1422)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.96165e-06,'m^3/(mol*s)'), n=3.30609, Ea=(150.918,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s',), comment="""Estimated from node CY_2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Br1s_5Cl1sF1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['C=C([O])C([O])C(F)F(4687)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=[C]C[CH]C(F)F(833)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=O(129)', '[CH2]C([O])C(F)F(4737)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['O=C1CC(C(F)F)O1(4695)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['O=CCC(=O)C(F)F(4835)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['O=C=CC(O)C(F)F(4836)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CO(13)', '[CH2]C([O])C(F)F(4737)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(27.525,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=O(338)', 'O=CC(F)F(990)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(46.8293,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHF2(82)', 'O=[C]CC=O(1325)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(37.3735,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'O=[C]CC(=O)C(F)F(4837)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(16.9476,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH2CO(28)', '[O][CH]C(F)F(1636)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(10.2251,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[O]C(C=C=O)C(F)F(4838)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00594739,'m^3/(mol*s)'), n=2.86123, Ea=(1.31722,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2779200861886376, var=1.6522143348846914, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=O(338)', '[O][CH]C(F)F(1636)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'O=[C]CC(=O)[CH]F(4839)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(245.963,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['O=[C]C[C](O)C(F)F(4840)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['[O]C=CC([O])C(F)F(4841)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(C[C]=O)C(F)F(4691)'],
    products = ['O=[C][CH]C(O)C(F)F(4842)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=[C]CC(O)[C](F)F(4843)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.248e+06,'s^-1'), n=0.992, Ea=(28.4919,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 441 used for R3H_SS_Cs;C_rad_out_noH;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_noH;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 27.7 to 28.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][C](CC=O)C(F)F(4844)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.98524e+06,'s^-1'), n=1.75122, Ea=(105.056,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(CC=O)[C](F)F(4845)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(446.55,'s^-1'), n=2.5049, Ea=(73.5615,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;XH_out] for rate rule [R4H_SSS;C_rad_out_noH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CC([CH]F)OF(4846)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(98.7168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C([CH]F)CC(=O)F(4847)'],
    products = ['[O]C(C[C]=O)C(F)F(4691)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(240.983,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

network(
    label = 'PDepNetwork #1890',
    isomers = [
        '[O]C(C[C]=O)C(F)F(4691)',
    ],
    reactants = [
        ('CH2CO(28)', 'O=CC(F)F(990)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1890',
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

