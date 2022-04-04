species(
    label = 'FC#CC(F)F(3950)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5 C u0 p0 c0 {4,S} {6,T}
6 C u0 p0 c0 {3,S} {5,T}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-325.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([143,454,630,930,1107,1168,1393,3045,2175,525,239,401,1367,347.84],'cm^-1')),
        HinderedRotor(inertia=(0.00139456,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.80767,0.015724,5.41417e-05,-1.64795e-07,1.29627e-10,-39107.5,10.8692], Tmin=(10,'K'), Tmax=(455.369,'K')), NASAPolynomial(coeffs=[4.56328,0.0230449,-1.59525e-05,5.13748e-09,-6.22916e-13,-39321,6.22913], Tmin=(455.369,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-325.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC#CC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C2HF(58)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[C]#CF(235)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u1 p0 c0 {2,T}
"""),
    E0 = (432.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([239,401,1367,1188.34],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0197,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2576.58,'J/mol'), sigma=(4.40467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=402.46 K, Pc=68.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81645,0.0171381,-3.27593e-05,2.86044e-08,-9.00898e-12,52027.9,8.68892], Tmin=(100,'K'), Tmax=(974.887,'K')), NASAPolynomial(coeffs=[4.85511,0.00251795,-6.39199e-07,5.75084e-11,-5.61704e-16,51927.7,0.430007], Tmin=(974.887,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.396,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Acetyl)"""),
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
    label = 'F[C]C#CF(10567)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,T}
4 C u0 p1 c0 {1,S} {3,S}
5 C u0 p0 c0 {2,S} {3,T}
"""),
    E0 = (235.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,252,948,239,401,1367,381.92],'cm^-1')),
        HinderedRotor(inertia=(0.00115494,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91232,0.0258665,-3.48928e-05,2.5621e-08,-7.66626e-12,28329,12.6533], Tmin=(100,'K'), Tmax=(811.239,'K')), NASAPolynomial(coeffs=[6.20239,0.00964357,-4.89515e-06,9.68485e-10,-6.88447e-14,27795.2,-2.53107], Tmin=(811.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtCs) + group(CJ2_singlet-FC) + group(CtCF)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'C2F2(60)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
"""),
    E0 = (-6.58462,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(61.9968,'amu')),
        LinearRotor(inertia=(141.582,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([284.265,284.265,295.939,295.939,811.153,1396.68,2594.35],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.0181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37365,0.00967198,6.42941e-05,-2.48371e-07,2.43362e-10,-791.12,6.86461], Tmin=(10,'K'), Tmax=(392.056,'K')), NASAPolynomial(coeffs=[5.56744,0.00756355,-5.20729e-06,1.71201e-09,-2.14827e-13,-1118.95,-3.65215], Tmin=(392.056,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-6.58462,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""FC#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC#C[CH]F(4859)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u1 p0 c0 {1,S} {4,S} {6,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u0 p0 c0 {2,S} {4,T}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (45.0609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([300,424,1135,1289,3214,2175,525,239,401,1367,3931.13],'cm^-1')),
        HinderedRotor(inertia=(0.0154395,'amu*angstrom^2'), symmetry=1, barrier=(18.3074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90024,0.00741053,7.60747e-05,-2.30944e-07,2.00261e-10,5421.02,10.0255], Tmin=(10,'K'), Tmax=(407.55,'K')), NASAPolynomial(coeffs=[4.65259,0.0154555,-1.03221e-05,3.27509e-09,-3.95979e-13,5231.56,5.49913], Tmin=(407.55,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(45.0609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), label="""FC#C[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[C]#CC(F)F(12621)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u1 p0 c0 {4,T}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (135.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([143,454,630,930,1107,1168,1393,3045,2175,525,188.573],'cm^-1')),
        HinderedRotor(inertia=(0.00468865,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.874,0.0191158,-1.299e-05,3.97785e-09,-4.11821e-13,16288.6,10.1534], Tmin=(10,'K'), Tmax=(1424,'K')), NASAPolynomial(coeffs=[8.99144,0.00730449,-3.24863e-06,6.81485e-10,-5.50496e-14,14571.3,-17.2568], Tmin=(1424,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(135.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""[C]#CC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]=C=C(F)F(4215)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {4,D} {6,D}
6 C u1 p0 c0 {3,S} {5,D}
"""),
    E0 = (-151.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([94,120,354,641,825,1294,540,610,2055,137,207,812],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81697,0.0142144,8.17672e-05,-2.81639e-07,2.55173e-10,-18237.7,11.4722], Tmin=(10,'K'), Tmax=(406.163,'K')), NASAPolynomial(coeffs=[5.51158,0.0188349,-1.39947e-05,4.71475e-09,-5.89764e-13,-18551.1,2.65979], Tmin=(406.163,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-151.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""F[C]DCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[C]=C(F)C(F)F-2(15445)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5 C u0 p0 c0 {3,S} {4,S} {6,D}
6 C u0 p1 c0 {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-215.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([171,205,598,1104,1143,1317,1411,3153,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.57959,'amu*angstrom^2'), symmetry=1, barrier=(36.318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85497,0.0521447,-8.85575e-05,7.8609e-08,-2.74463e-11,-25850.4,15.8141], Tmin=(100,'K'), Tmax=(795.168,'K')), NASAPolynomial(coeffs=[7.71183,0.0163119,-8.94557e-06,1.78728e-09,-1.26002e-13,-26580.4,-9.83312], Tmin=(795.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCCF) + group(CdJ2_singlet-Cds)"""),
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
    E0 = (5.5112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-55.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (267.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (110.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (197.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (165.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (51.8577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-58.4945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CF2(43)', 'C2HF(58)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.44768e+28,'m^3/(mol*s)'), n=-6.4458, Ea=(124.34,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.6591608693425184, var=5.4995193120720405, Tref=1000.0, N=19, data_mean=0.0, correlation='CH',), comment="""Estimated from node CH"""),
)

reaction(
    label = 'reaction2',
    reactants = ['HF(38)', 'F[C]C#CF(10567)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.9697e-08,'m^3/(mol*s)'), n=3.78805, Ea=(0.793798,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3001410487914693, var=2.9282398409420045, Tref=1000.0, N=3, data_mean=0.0, correlation='HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R',), comment="""Estimated from node HY_N-3Br1sCCl1sF1sHI1s->F1s_N-4Br1sCl1sF1sI1s->Br1s_N-2Br1sCl1sF1sHI1s->H_Ext-3Br1sCCl1sHI1s-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'C2F2(60)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.00938e+53,'m^3/(mol*s)'), n=-13.541, Ea=(145.965,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'FC#C[CH]F(4859)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(3.01345,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', '[C]#CC(F)F(12621)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.38619e+06,'m^3/(mol*s)'), n=0.213797, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.00016139601674549603, var=0.035561666158317407, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-3BrCO-R_N-Sp-3BrBrCCOO=1C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHF2(82)', '[C]#CF(235)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -42.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'F[C]=C=C(F)F(4215)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(2.13853,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]=C(F)C(F)F-2(15445)'],
    products = ['FC#CC(F)F(3950)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(167.488,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

network(
    label = 'PDepNetwork #870',
    isomers = [
        'FC#CC(F)F(3950)',
    ],
    reactants = [
        ('CF2(43)', 'C2HF(58)'),
        ('CHF2(82)', '[C]#CF(235)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #870',
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

