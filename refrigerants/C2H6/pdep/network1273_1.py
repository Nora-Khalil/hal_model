species(
    label = 'C=CC(=C)C(=O)O(1876)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {3,S} {7,D} {8,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {3,D} {9,S} {10,S}
7  C u0 p0 c0 {4,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-192.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,381.172,381.252,381.261,381.632,381.682,381.842],'cm^-1')),
        HinderedRotor(inertia=(0.00115603,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102548,'amu*angstrom^2'), symmetry=1, barrier=(10.5978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102485,'amu*angstrom^2'), symmetry=1, barrier=(10.6008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922231,0.0619498,-5.7623e-05,2.72634e-08,-5.11201e-12,-23041.4,23.8184], Tmin=(100,'K'), Tmax=(1290.24,'K')), NASAPolynomial(coeffs=[14.83,0.0188329,-7.49656e-06,1.36313e-09,-9.35264e-14,-26630.3,-46.823], Tmin=(1290.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CO2(15)',
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
    label = 'C=C=CC(1692)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.759585,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74634,0.0218191,8.22302e-06,-2.14762e-08,8.55596e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.61,'K')), NASAPolynomial(coeffs=[6.82084,0.0192337,-7.45616e-06,1.36535e-09,-9.53184e-14,16028,-10.4336], Tmin=(1025.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC=C(61)',
    structure = adjacencyList("""1  C u0 p0 c0 {2,S} {3,D} {5,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {7,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80601,0.0102582,6.17268e-05,-9.01653e-08,3.59121e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.045,'K')), NASAPolynomial(coeffs=[12.4693,0.0100555,-2.41212e-06,4.5709e-10,-3.93171e-14,8010.8,-43.6372], Tmin=(946.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CO(14)',
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
    label = 'C=CC(=C)O(570)',
    structure = adjacencyList("""1  O u0 p2 c0 {2,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {6,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {7,S} {8,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {1,S}
"""),
    E0 = (-119.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.0349,'amu*angstrom^2'), symmetry=1, barrier=(23.7943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04062,'amu*angstrom^2'), symmetry=1, barrier=(23.9259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62495,0.0398226,-1.28515e-06,-3.80372e-08,2.09279e-11,-14332.2,15.7636], Tmin=(100,'K'), Tmax=(917.013,'K')), NASAPolynomial(coeffs=[16.745,0.00649243,-1.28844e-07,-8.27573e-11,4.0961e-15,-18476.9,-63.3508], Tmin=(917.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=C=O)CO(2723)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  C u0 p0 c0 {2,D} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-162.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.420299,'amu*angstrom^2'), symmetry=1, barrier=(9.66351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.420394,'amu*angstrom^2'), symmetry=1, barrier=(9.66569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.420527,'amu*angstrom^2'), symmetry=1, barrier=(9.66875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890028,0.0762609,-0.000118601,1.04888e-07,-3.62566e-11,-19443.6,24.0444], Tmin=(100,'K'), Tmax=(858.467,'K')), NASAPolynomial(coeffs=[6.69661,0.0331518,-1.52264e-05,2.82655e-09,-1.90899e-13,-19849,0.362426], Tmin=(858.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cd-Cd(CCO)H) + group(Cds-CdsHH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'C=C=CCC(=O)O(1886)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,D} {10,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u0 p0 c0 {4,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-218.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,180,180,606.489,888.84,889.29,889.743],'cm^-1')),
        HinderedRotor(inertia=(1.15698,'amu*angstrom^2'), symmetry=1, barrier=(26.6012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473868,'amu*angstrom^2'), symmetry=1, barrier=(26.601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15696,'amu*angstrom^2'), symmetry=1, barrier=(26.6008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4269.97,'J/mol'), sigma=(6.32747,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=666.96 K, Pc=38.25 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46076,0.046517,-1.52808e-05,-8.58872e-09,4.68425e-12,-26150,22.9668], Tmin=(100,'K'), Tmax=(1216.95,'K')), NASAPolynomial(coeffs=[13.4683,0.0256834,-1.25696e-05,2.50839e-09,-1.80261e-13,-30452.3,-42.9895], Tmin=(1216.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=C(O)OC=C(2724)',
    structure = adjacencyList("""1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {3,S} {13,S}
3  C u0 p0 c0 {1,S} {2,S} {7,D}
4  C u0 p0 c0 {1,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-46.7442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.965697,'amu*angstrom^2'), symmetry=1, barrier=(22.2033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.967929,'amu*angstrom^2'), symmetry=1, barrier=(22.2546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964324,'amu*angstrom^2'), symmetry=1, barrier=(22.1717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440924,0.0715178,-7.26521e-05,3.47397e-08,-5.79024e-12,-5487.94,24.4598], Tmin=(100,'K'), Tmax=(1006.9,'K')), NASAPolynomial(coeffs=[17.3047,0.0154258,-5.33005e-06,9.18125e-10,-6.24264e-14,-9436.57,-59.7583], Tmin=(1006.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.7442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C1=C(O)OC[CH]1(2725)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u1 p0 c0 {3,S} {4,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {4,D}
7  C u1 p0 c0 {4,S} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
"""),
    E0 = (-56.2996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905882,0.0490685,1.50744e-05,-7.34616e-08,3.78887e-11,-6642.01,15.4709], Tmin=(100,'K'), Tmax=(899.862,'K')), NASAPolynomial(coeffs=[22.3183,0.005256,2.47914e-06,-6.9302e-10,4.79375e-14,-12575.4,-97.1286], Tmin=(899.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.2996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(2,3-Dihydrofuran) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = 'O=C(O)C1=CCC1(2320)',
    structure = adjacencyList("""1  O u0 p2 c0 {7,S} {13,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {12,S}
7  C u0 p0 c0 {1,S} {2,D} {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (-158.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7926,0.0412294,-1.34955e-05,-6.09168e-09,3.5463e-12,-19019.1,23.3957], Tmin=(100,'K'), Tmax=(1196.98,'K')), NASAPolynomial(coeffs=[10.4392,0.0258971,-1.12777e-05,2.139e-09,-1.49781e-13,-22060.7,-23.9325], Tmin=(1196.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cyclobutene)"""),
)

species(
    label = 'C=CC(=C)C([O])[O](2726)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (147.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,834.044],'cm^-1')),
        HinderedRotor(inertia=(0.121978,'amu*angstrom^2'), symmetry=1, barrier=(2.80452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120918,'amu*angstrom^2'), symmetry=1, barrier=(2.78015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3485,0.0639279,-8.28456e-05,7.12735e-08,-2.62468e-11,17852.2,25.0051], Tmin=(100,'K'), Tmax=(749.638,'K')), NASAPolynomial(coeffs=[5.25088,0.0373591,-1.81847e-05,3.54433e-09,-2.49524e-13,17428.6,8.37969], Tmin=(749.638,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(=O)O(2035)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,D} {3,S}
5  C u1 p0 c0 {3,S} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (52.6694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00596,0.0716709,-9.98712e-05,8.4304e-08,-2.93374e-11,6436.9,26.4289], Tmin=(100,'K'), Tmax=(784.631,'K')), NASAPolynomial(coeffs=[7.09913,0.0343082,-1.63998e-05,3.14868e-09,-2.19012e-13,5674.65,-0.253395], Tmin=(784.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.6694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(=C)C([O])O(2727)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {13,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {7,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (120.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,184.636,184.733,184.794],'cm^-1')),
        HinderedRotor(inertia=(0.00493694,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.512881,'amu*angstrom^2'), symmetry=1, barrier=(12.4055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.51863,'amu*angstrom^2'), symmetry=1, barrier=(85.203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85455,0.0738193,-0.000101462,8.066e-08,-2.60464e-11,14658.5,26.2083], Tmin=(100,'K'), Tmax=(813.995,'K')), NASAPolynomial(coeffs=[8.92638,0.0301146,-1.3481e-05,2.50652e-09,-1.70981e-13,13478.2,-10.2504], Tmin=(813.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(C[CH2])C(=O)O(2728)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {1,S} {2,D} {4,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u1 p0 c0 {4,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (132.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07193,0.0527493,-1.46899e-05,-9.80741e-08,1.1264e-10,15995.7,25.565], Tmin=(100,'K'), Tmax=(442.758,'K')), NASAPolynomial(coeffs=[6.29597,0.0339241,-1.64211e-05,3.16884e-09,-2.20739e-13,15432.1,6.48764], Tmin=(442.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C([O])O(2729)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {12,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {9,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  C u1 p0 c0 {4,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (169.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,339.276,339.32],'cm^-1')),
        HinderedRotor(inertia=(0.141393,'amu*angstrom^2'), symmetry=1, barrier=(11.5469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141381,'amu*angstrom^2'), symmetry=1, barrier=(11.5464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141149,'amu*angstrom^2'), symmetry=1, barrier=(11.5463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07194,0.0681248,-7.88794e-05,5.07527e-08,-1.34257e-11,20437.1,26.2452], Tmin=(100,'K'), Tmax=(910.502,'K')), NASAPolynomial(coeffs=[10.1783,0.0281189,-1.29722e-05,2.49575e-09,-1.75653e-13,18778.8,-16.8341], Tmin=(910.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(=O)O(2037)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  C u1 p0 c0 {4,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (61.9238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09408,0.0688836,-8.47633e-05,6.17088e-08,-1.88553e-11,7547.83,26.054], Tmin=(100,'K'), Tmax=(787.517,'K')), NASAPolynomial(coeffs=[8.26526,0.0324585,-1.53818e-05,2.97294e-09,-2.08961e-13,6418.38,-6.82964], Tmin=(787.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)C([O])=O(1879)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u0 p0 c0 {4,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (40.5328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4196.76,'J/mol'), sigma=(6.40471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=655.52 K, Pc=36.25 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38219,0.064598,-8.86303e-05,8.23345e-08,-3.17795e-11,4962.41,24.7694], Tmin=(100,'K'), Tmax=(782.472,'K')), NASAPolynomial(coeffs=[3.4969,0.0413989,-2.04086e-05,3.97528e-09,-2.78846e-13,5010.72,17.5093], Tmin=(782.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.5328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([C]=C)=C(O)O(2730)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {4,S} {12,S}
3  C u0 p0 c0 {4,D} {5,S} {7,S}
4  C u0 p0 c0 {1,S} {2,S} {3,D}
5  C u1 p0 c0 {3,S} {8,S} {9,S}
6  C u0 p0 c0 {7,D} {10,S} {11,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (40.8959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.18425,'amu*angstrom^2'), symmetry=1, barrier=(27.2282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18407,'amu*angstrom^2'), symmetry=1, barrier=(27.224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18438,'amu*angstrom^2'), symmetry=1, barrier=(27.2313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18394,'amu*angstrom^2'), symmetry=1, barrier=(27.221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.789092,0.0899362,-0.000108334,6.10559e-08,-1.26972e-11,5104.37,25.2996], Tmin=(100,'K'), Tmax=(1358.11,'K')), NASAPolynomial(coeffs=[23.5994,0.00413139,1.86879e-06,-6.16564e-10,5.00345e-14,-231.362,-95.0816], Tmin=(1358.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.8959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C(=CC)C(=O)O(2670)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {4,D} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  C u2 p0 c0 {4,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (37.4709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29285,0.0659319,-8.35077e-05,7.43498e-08,-2.81651e-11,4598,25.9571], Tmin=(100,'K'), Tmax=(783.955,'K')), NASAPolynomial(coeffs=[3.57943,0.0437466,-2.09333e-05,4.02257e-09,-2.80276e-13,4562.7,17.5438], Tmin=(783.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.4709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)O2s) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C)=C(O)O(2731)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,S} {11,S}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {3,D}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u2 p0 c0 {3,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (61.0857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.94355,'amu*angstrom^2'), symmetry=1, barrier=(44.686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94389,'amu*angstrom^2'), symmetry=1, barrier=(44.6939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94376,'amu*angstrom^2'), symmetry=1, barrier=(44.6908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94183,'amu*angstrom^2'), symmetry=1, barrier=(44.6466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.896353,0.0878154,-9.35819e-05,4.78145e-08,-9.11582e-12,7540.24,26.9283], Tmin=(100,'K'), Tmax=(1475.53,'K')), NASAPolynomial(coeffs=[23.2689,0.00934865,-6.41729e-07,-1.2865e-10,1.55503e-14,1819.46,-94.2761], Tmin=(1475.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.0857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[C](C)C(=O)O(2732)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {6,S} {13,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  C u1 p0 c0 {5,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (3.98694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72966,0.0515528,-3.62468e-05,1.31681e-08,-2.03176e-12,559.413,23.5972], Tmin=(100,'K'), Tmax=(1432.31,'K')), NASAPolynomial(coeffs=[9.54129,0.0297373,-1.34002e-05,2.53426e-09,-1.75686e-13,-1678.32,-16.896], Tmin=(1432.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.98694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCJ(C)CO) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)=C([O])[O](2033)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {1,S} {2,S} {4,D}
7  C u0 p0 c0 {5,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-26.6735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.879681,'amu*angstrom^2'), symmetry=1, barrier=(20.2256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.878178,'amu*angstrom^2'), symmetry=1, barrier=(20.191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675511,0.065171,-5.66375e-05,2.1387e-08,-2.05914e-12,-3081.48,22.9953], Tmin=(100,'K'), Tmax=(1035.58,'K')), NASAPolynomial(coeffs=[16.2869,0.0177512,-6.60793e-06,1.19011e-09,-8.27912e-14,-7005.5,-56.201], Tmin=(1035.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.6735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC(=C)C([O])O(2733)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {3,S} {12,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u1 p0 c0 {5,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (169.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,339.276,339.32],'cm^-1')),
        HinderedRotor(inertia=(0.141393,'amu*angstrom^2'), symmetry=1, barrier=(11.5469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141381,'amu*angstrom^2'), symmetry=1, barrier=(11.5464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141149,'amu*angstrom^2'), symmetry=1, barrier=(11.5463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07194,0.0681248,-7.88794e-05,5.07527e-08,-1.34257e-11,20437.1,26.2452], Tmin=(100,'K'), Tmax=(910.502,'K')), NASAPolynomial(coeffs=[10.1783,0.0281189,-1.29722e-05,2.49575e-09,-1.75653e-13,18778.8,-16.8341], Tmin=(910.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC(=C)C([O])=O(2668)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {7,D}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  C u0 p0 c0 {4,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (111.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,202.728,807.926,968.735,1128.35,1303.85,1492.41,1653.25],'cm^-1')),
        HinderedRotor(inertia=(0.149248,'amu*angstrom^2'), symmetry=1, barrier=(3.49423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149248,'amu*angstrom^2'), symmetry=1, barrier=(3.49423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149248,'amu*angstrom^2'), symmetry=1, barrier=(3.49423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81217,0.0376244,-1.72424e-05,2.73656e-09,-1.13481e-13,13336,18.6687], Tmin=(100,'K'), Tmax=(2821.19,'K')), NASAPolynomial(coeffs=[44.8255,-0.00886468,1.27523e-06,-1.74118e-10,1.4613e-14,-14445.9,-229.959], Tmin=(2821.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC([CH2])=C(O)O(2734)',
    structure = adjacencyList("""multiplicity 3
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {4,S} {11,S}
3  C u0 p0 c0 {4,D} {5,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {3,D}
5  C u0 p0 c0 {3,S} {7,D} {8,S}
6  C u1 p0 c0 {3,S} {9,S} {10,S}
7  C u1 p0 c0 {5,D} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (88.9966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.05308,'amu*angstrom^2'), symmetry=1, barrier=(24.2123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04582,'amu*angstrom^2'), symmetry=1, barrier=(24.0454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04569,'amu*angstrom^2'), symmetry=1, barrier=(24.0424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04716,'amu*angstrom^2'), symmetry=1, barrier=(24.0764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00978,0.089457,-0.000104291,5.5774e-08,-1.09207e-11,10902,26.9054], Tmin=(100,'K'), Tmax=(1460.69,'K')), NASAPolynomial(coeffs=[25.2926,0.00122007,2.96631e-06,-7.75988e-10,5.82515e-14,4947.34,-104.035], Tmin=(1460.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.9966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=CC)C([O])=O(2445)',
    structure = adjacencyList("""multiplicity 3
1  O u1 p2 c0 {7,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,D} {6,S} {7,S}
5  C u0 p0 c0 {3,S} {4,D} {11,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {2,D} {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (48.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,180,180,1684.59,1685.24,1685.26],'cm^-1')),
        HinderedRotor(inertia=(0.149161,'amu*angstrom^2'), symmetry=1, barrier=(3.4295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148262,'amu*angstrom^2'), symmetry=1, barrier=(3.40884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148119,'amu*angstrom^2'), symmetry=1, barrier=(3.40555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4276.76,'J/mol'), sigma=(6.46324,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=668.02 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53531,0.062329,-9.08789e-05,8.87556e-08,-3.49488e-11,5955.29,23.6425], Tmin=(100,'K'), Tmax=(800.839,'K')), NASAPolynomial(coeffs=[2.26886,0.0417762,-2.07492e-05,4.04175e-09,-2.82824e-13,6379.38,23.6478], Tmin=(800.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)O2s) + radical(CCOJ) + radical(C=C(C=O)CJ)"""),
)

species(
    label = 'OH(8)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39807e-08,-2.13441e-11,2.48061e-15,3579.39,4.57802], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=CC(=C)[C]=O(2365)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,D} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (143.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.07471,'amu*angstrom^2'), symmetry=1, barrier=(24.7096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07027,'amu*angstrom^2'), symmetry=1, barrier=(24.6077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0925,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32704,0.0495299,-3.19523e-05,9.86414e-10,4.13383e-12,17372,19.3902], Tmin=(100,'K'), Tmax=(1030.42,'K')), NASAPolynomial(coeffs=[15.5037,0.0128629,-5.31065e-06,1.04709e-09,-7.76085e-14,13475.4,-54.1597], Tmin=(1030.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'H(6)',
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
    label = 'C=CC(=C)C([O])=O(2032)',
    structure = adjacencyList("""multiplicity 2
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,D} {6,S}
4  C u0 p0 c0 {3,S} {7,D} {8,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {2,D} {3,S}
7  C u0 p0 c0 {4,D} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (33.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,273.69,273.691,273.691,273.691,273.696,757.499,4000],'cm^-1')),
        HinderedRotor(inertia=(0.367231,'amu*angstrom^2'), symmetry=1, barrier=(19.5197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126713,'amu*angstrom^2'), symmetry=1, barrier=(19.5197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0919,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72115,0.0539357,-5.09978e-05,2.71786e-08,-6.21273e-12,4067.96,21.1276], Tmin=(100,'K'), Tmax=(1017.79,'K')), NASAPolynomial(coeffs=[8.18499,0.0285322,-1.35586e-05,2.65529e-09,-1.89061e-13,2752.19,-10.1706], Tmin=(1017.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'C2H3(28)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=[C]C(=O)O(2735)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {8,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {5,S}
4 C u0 p0 c0 {5,D} {6,S} {7,S}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-43.6283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1685,370,180,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.259459,'amu*angstrom^2'), symmetry=1, barrier=(5.96548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00052067,'amu*angstrom^2'), symmetry=1, barrier=(5.91172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98535,0.0251532,-3.97798e-05,3.82837e-08,-1.50862e-11,-5213.49,14.5134], Tmin=(100,'K'), Tmax=(727.078,'K')), NASAPolynomial(coeffs=[4.59246,0.0127631,-6.89751e-06,1.42079e-09,-1.03072e-13,-5353.39,7.91732], Tmin=(727.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.6283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=CJC=O)"""),
)

species(
    label = 'O=[C]O(174)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (-192.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1244.76,1684.94],'cm^-1')),
        HinderedRotor(inertia=(0.00103939,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-192.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=[C]C=C(1095)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {4,D} {8,S} {9,S}
4 C u1 p0 c0 {1,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (303.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.39164,'amu*angstrom^2'), symmetry=1, barrier=(31.9966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0824,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67736,0.0217333,1.0387e-05,-2.87121e-08,1.2538e-11,36547.1,13.0516], Tmin=(100,'K'), Tmax=(966.654,'K')), NASAPolynomial(coeffs=[8.97812,0.0135587,-4.70087e-06,8.4736e-10,-6.05066e-14,34492.7,-21.4574], Tmin=(966.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(=C=C)C(=O)O(2736)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {12,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {5,S} {7,D}
4  C u0 p0 c0 {1,S} {2,D} {3,S}
5  C u1 p0 c0 {3,S} {8,S} {9,S}
6  C u0 p0 c0 {7,D} {10,S} {11,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
"""),
    E0 = (-2.35716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,345.801,346.265,346.55,347.261,347.549],'cm^-1')),
        HinderedRotor(inertia=(0.00141051,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140649,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160898,'amu*angstrom^2'), symmetry=1, barrier=(13.6541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0919,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19747,0.063064,-6.93823e-05,4.03021e-08,-9.44584e-12,-183.715,22.5233], Tmin=(100,'K'), Tmax=(1030.38,'K')), NASAPolynomial(coeffs=[11.6912,0.0223268,-1.00782e-05,1.93194e-09,-1.36159e-13,-2346.23,-28.4171], Tmin=(1030.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.35716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH]=C(C=C)C(=O)O(2737)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {5,S} {11,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {7,D}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {4,D} {9,S} {10,S}
7  C u1 p0 c0 {3,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (54.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,344.387,344.389,344.391,344.391,344.394],'cm^-1')),
        HinderedRotor(inertia=(0.00142135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1325,'amu*angstrom^2'), symmetry=1, barrier=(11.1514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132497,'amu*angstrom^2'), symmetry=1, barrier=(11.1514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0919,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00169,0.0636925,-6.82214e-05,3.70915e-08,-7.95739e-12,6671.73,23.9353], Tmin=(100,'K'), Tmax=(1135.83,'K')), NASAPolynomial(coeffs=[14.0851,0.0176163,-7.37088e-06,1.37491e-09,-9.58742e-14,3699.7,-40.8505], Tmin=(1135.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(=C)C(=O)O(2738)',
    structure = adjacencyList("""multiplicity 2
1  O u0 p2 c0 {4,S} {11,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {5,S} {6,D}
4  C u0 p0 c0 {1,S} {2,D} {3,S}
5  C u0 p0 c0 {3,S} {7,D} {8,S}
6  C u0 p0 c0 {3,D} {9,S} {10,S}
7  C u1 p0 c0 {5,D} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (54.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,344.387,344.389,344.391,344.391,344.394],'cm^-1')),
        HinderedRotor(inertia=(0.00142135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1325,'amu*angstrom^2'), symmetry=1, barrier=(11.1514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132497,'amu*angstrom^2'), symmetry=1, barrier=(11.1514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0919,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00169,0.0636925,-6.82214e-05,3.70915e-08,-7.95739e-12,6671.73,23.9353], Tmin=(100,'K'), Tmax=(1135.83,'K')), NASAPolynomial(coeffs=[14.0851,0.0176163,-7.37088e-06,1.37491e-09,-9.58742e-14,3699.7,-40.8505], Tmin=(1135.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(22)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C=C(O)O(2739)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {8,S}
2 O u0 p2 c0 {3,S} {9,S}
3 C u0 p0 c0 {1,S} {2,S} {5,D}
4 C u0 p0 c0 {5,D} {6,S} {7,S}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-183.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0626,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889118,0.0557305,-6.5491e-05,3.53566e-08,-6.97911e-12,-21911.1,16.3357], Tmin=(100,'K'), Tmax=(1452.9,'K')), NASAPolynomial(coeffs=[17.1073,0.000717385,2.00405e-06,-5.22732e-10,3.93408e-14,-25530.1,-64.2024], Tmin=(1452.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C([C]C)C(=O)O(2740)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u0 p0 c0 {1,S} {2,D} {4,S}
6  C u0 p0 c0 {4,D} {11,S} {12,S}
7  C u0 p1 c0 {3,S} {4,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (111.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966525,0.0770226,-0.000125306,1.19603e-07,-4.4454e-11,13529.7,35.2923], Tmin=(100,'K'), Tmax=(833.542,'K')), NASAPolynomial(coeffs=[4.00063,0.0405122,-2.01037e-05,3.86893e-09,-2.67189e-13,13786.4,25.781], Tmin=(833.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(C=C)C(=O)O-2(2741)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {1,S} {2,D} {3,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p1 c0 {3,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (53.9261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76696,0.052007,-3.72843e-05,1.37088e-08,-2.14578e-12,6563.28,23.7423], Tmin=(100,'K'), Tmax=(1406.24,'K')), NASAPolynomial(coeffs=[9.4187,0.0302418,-1.40679e-05,2.70242e-09,-1.89064e-13,4411.25,-15.7815], Tmin=(1406.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.9261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC(=C)C(=O)O(2742)',
    structure = adjacencyList("""1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {1,S} {2,D} {4,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  C u0 p1 c0 {3,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
"""),
    E0 = (130.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180,180,980.219,1600,1967.14,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.107006,'amu*angstrom^2'), symmetry=1, barrier=(2.46028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107006,'amu*angstrom^2'), symmetry=1, barrier=(2.46028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107006,'amu*angstrom^2'), symmetry=1, barrier=(2.46028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107006,'amu*angstrom^2'), symmetry=1, barrier=(2.46028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0998,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11245,0.073687,-0.000121704,1.17797e-07,-4.43639e-11,15745.3,26.0725], Tmin=(100,'K'), Tmax=(825.553,'K')), NASAPolynomial(coeffs=[3.77448,0.039203,-1.98269e-05,3.85491e-09,-2.67945e-13,16041.4,18.1951], Tmin=(825.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CsCsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-113.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-7.88242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-22.1998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (19.7499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-77.5476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (17.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-32.2759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-102.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (138.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (43.0754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (177.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (163.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (200.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (92.8685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (71.4775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (69.1085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (29.9887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (63.1483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-10.6865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-41.3469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (161.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (103.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (91.0592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (41.3609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (139.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (212.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (210.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (78.8699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (177.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (233.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (233.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (57.3176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (106.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (48.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (111.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=CC(=C)C(=O)O(1876)'],
    products = ['CO2(15)', 'C=C=CC(1692)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(111.899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(14)', 'C=CC(=C)O(570)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(263.285,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CO2(15)', 'C=CC=C(61)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(20.4811,'m^3/(mol*s)'), n=1.86833, Ea=(316.938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO2_Cdd;R_H] for rate rule [CO2_Cdd;Cd/H/OneDe]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: 1,3_Insertion_CO2"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC(=C=O)CO(2723)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(214.737,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CC(=C)C(=O)O(1876)'],
    products = ['C=C=CCC(=O)O(1886)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(147.446,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C=C(O)OC=C(2724)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(97.0146,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C1=C(O)OC[CH]1(2725)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.43612e+09,'s^-1'), n=0.0758668, Ea=(56.4791,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC(=C)C(=O)O(1876)'],
    products = ['O=C(O)C1=CCC1(2320)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH2_1;CdH2_2]
Euclidian distance = 1.4142135623730951
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC(=C)C([O])[O](2726)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([C]=C)C(=O)O(2035)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(=C)C([O])O(2727)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C[CH2])C(=O)O(2728)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=C)C([O])O(2729)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC([CH2])C(=O)O(2037)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(C=C)C([O])=O(1879)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([C]=C)=C(O)O(2730)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(=CC)C(=O)O(2670)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(C=C)=C(O)O(2731)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C[C](C)C(=O)O(2732)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.77e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=CC(C)=C([O])[O](2033)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.554e+09,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC(=C)C([O])O(2733)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CC(=C)C([O])=O(2668)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=CC([CH2])=C(O)O(2734)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=CC)C([O])=O(2445)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.27566e+10,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['OH(8)', 'C=CC(=C)[C]=O(2365)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R
Ea raised from -1.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(6)', 'C=CC(=C)C([O])=O(2032)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C2H3(28)', 'C=[C]C(=O)O(2735)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=[C]O(174)', 'C=[C]C=C(1095)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(6)', '[CH2]C(=C=C)C(=O)O(2736)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(0.143818,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(6)', '[CH]=C(C=C)C(=O)O(2737)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -5.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(6)', '[CH]=CC(=C)C(=O)O(2738)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.83712e+08,'m^3/(mol*s)'), n=-0.251115, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.04900385874705703, var=3.2086570313419123, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_4R!H->C_N-Sp-3C-2C
Ea raised from -5.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=CC(=C)C(=O)O(1876)'],
    products = ['C2H2(22)', 'C=C=C(O)O(2739)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(23.1678,'s^-1'), n=3.09547, Ea=(282.311,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7BrCClFIOPSSi->O',), comment="""Estimated from node Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7BrCClFIOPSSi->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([C]C)C(=O)O(2740)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.99999e+12,'s^-1'), n=8.2394e-08, Ea=(27.5967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_N-1C-inRing
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C(C=C)C(=O)O-2(2741)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(27.3103,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]CC(=C)C(=O)O(2742)'],
    products = ['C=CC(=C)C(=O)O(1876)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.40886e+17,'s^-1'), n=-1.43042, Ea=(13.9361,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1273',
    isomers = [
        'C=CC(=C)C(=O)O(1876)',
    ],
    reactants = [
        ('CO2(15)', 'C=C=CC(1692)'),
        ('CO2(15)', 'C=CC=C(61)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1273',
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

