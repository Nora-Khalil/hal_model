species(
    label = 'O=CC(F)(F)OOOC(F)[CH]F(5994)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {12,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {11,S} {13,S}
11 C u1 p0 c0 {4,S} {10,S} {15,S}
12 C u0 p0 c0 {8,D} {9,S} {14,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-906.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.933142,0.115412,-0.000149387,9.68734e-08,-2.50021e-11,-108795,39.283], Tmin=(100,'K'), Tmax=(942.81,'K')), NASAPolynomial(coeffs=[18.6178,0.0324654,-1.74207e-05,3.55963e-09,-2.58798e-13,-112482,-53.8879], Tmin=(942.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-906.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = '[O]OC(F)(F)C=O(1651)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-551.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,251,367,519,700,855,1175,1303,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.588033,'amu*angstrom^2'), symmetry=1, barrier=(13.52,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.588634,'amu*angstrom^2'), symmetry=1, barrier=(13.5339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3731.44,'J/mol'), sigma=(6.05191,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.84 K, Pc=38.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5719,0.0584797,-9.64597e-05,8.34326e-08,-2.83599e-11,-66202.3,20.9945], Tmin=(100,'K'), Tmax=(811.123,'K')), NASAPolynomial(coeffs=[8.46089,0.0179335,-9.32175e-06,1.82202e-09,-1.26843e-13,-67103.6,-9.46575], Tmin=(811.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-551.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'FC1OC1F(977)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-489.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([205,317,423,563,514,686,1105,1199,1293,1437,1381,1463,3031,3163,1635.66],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2878.14,'J/mol'), sigma=(4.78691,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=449.56 K, Pc=59.54 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02952,-0.0033861,0.000100719,-1.79462e-07,1.00462e-10,-58841.1,9.38768], Tmin=(10,'K'), Tmax=(577.706,'K')), NASAPolynomial(coeffs=[2.18077,0.0257049,-1.71133e-05,5.32518e-09,-6.26769e-13,-58899.3,14.9396], Tmin=(577.706,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-489.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1OC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)[CH]F(975)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-236.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(1.00516,'amu*angstrom^2'), symmetry=1, barrier=(23.1107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.57,'J/mol'), sigma=(5.31812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.11 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18999,0.0346621,-3.10228e-05,1.06717e-08,-5.11859e-13,-28427.3,16.7377], Tmin=(100,'K'), Tmax=(1040.82,'K')), NASAPolynomial(coeffs=[12.0464,0.00571272,-2.17142e-06,4.35013e-10,-3.30333e-14,-30962.7,-33.5316], Tmin=(1040.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-236.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = 'F[CH]C(F)OOOC(F)F(5928)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {5,S} {6,S}
8  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {12,S}
10 C u1 p0 c0 {4,S} {8,S} {13,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-797.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,519,593,830,1115,1166,1389,1413,3114,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0635703,0.0911797,-0.000111594,6.86208e-08,-1.68414e-11,-95740.1,33.4842], Tmin=(100,'K'), Tmax=(988.245,'K')), NASAPolynomial(coeffs=[15.9293,0.0269628,-1.4124e-05,2.86918e-09,-2.08288e-13,-98876,-42.8715], Tmin=(988.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-797.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFHO) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'O=C[C]F-2(1630)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (35.6539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,262,1290],'cm^-1')),
        HinderedRotor(inertia=(0.407026,'amu*angstrom^2'), symmetry=1, barrier=(9.35834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79853,0.0176142,-2.7066e-05,3.1192e-08,-1.61007e-11,4288.74,9.0865], Tmin=(10,'K'), Tmax=(518.444,'K')), NASAPolynomial(coeffs=[4.38817,0.0119962,-7.71993e-06,2.33921e-09,-2.7033e-13,4241.97,6.76768], Tmin=(518.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(35.6539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(F)OOOF(5927)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {6,S}
6  O u0 p2 c0 {3,S} {5,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u1 p0 c0 {2,S} {7,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-255.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848142,0.0699678,-8.5007e-05,4.99213e-08,-1.14685e-11,-30558.7,27.2263], Tmin=(100,'K'), Tmax=(1064.19,'K')), NASAPolynomial(coeffs=[15.4067,0.0152453,-7.87305e-06,1.59964e-09,-1.16525e-13,-33657.3,-43.9158], Tmin=(1064.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2sFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)OOOC(F)=COF(6654)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {4,S} {12,S}
9  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {6,S} {12,D}
11 C u1 p0 c0 {3,S} {9,S} {15,S}
12 C u0 p0 c0 {8,S} {10,D} {14,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {12,S}
15 H u0 p0 c0 {11,S}
"""),
    E0 = (-454.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,487,638,688,1119,1325,1387,3149,326,540,652,719,1357,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07687,0.121932,-0.000180318,1.39459e-07,-4.33566e-11,-54481.2,41.8799], Tmin=(100,'K'), Tmax=(784.888,'K')), NASAPolynomial(coeffs=[15.2607,0.0386721,-2.12018e-05,4.31166e-09,-3.10369e-13,-57045.9,-32.9826], Tmin=(784.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-454.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2sCF) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[CH]C(F)OOOOC=C(F)F(6655)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {6,S} {7,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u1 p0 c0 {2,S} {9,S} {15,S}
11 C u0 p0 c0 {6,S} {12,D} {14,S}
12 C u0 p0 c0 {3,S} {4,S} {11,D}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-566.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.820337,0.105877,-0.000121642,6.79477e-08,-1.49106e-11,-67953.9,43.255], Tmin=(100,'K'), Tmax=(1109.76,'K')), NASAPolynomial(coeffs=[21.5086,0.0253942,-1.28575e-05,2.59674e-09,-1.88622e-13,-72909.8,-66.7943], Tmin=(1109.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-566.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = '[CH]C(F)OOOC(F)(F)C=O(6656)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
9  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
10 C u0 p0 c0 {7,D} {8,S} {13,S}
11 C u2 p0 c0 {9,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-482.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,261,493,600,1152,1365,1422,3097,2782.5,750,1395,475,1775,1000,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.404439,0.101094,-0.000120986,7.16739e-08,-1.68876e-11,-57848.3,36.721], Tmin=(100,'K'), Tmax=(1028.9,'K')), NASAPolynomial(coeffs=[18.2692,0.0284981,-1.5153e-05,3.10101e-09,-2.26173e-13,-61691.1,-53.9009], Tmin=(1028.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-482.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]F(223)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.000263365,8.89188e-06,-1.03032e-08,3.5081e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164132,-7.73117e-07,1.90988e-10,-1.59926e-14,25413.4,-0.815515], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'O=CC(F)(F)OOO[CH]F(3089)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {9,D}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
9  C u0 p0 c0 {7,D} {8,S} {11,S}
10 C u1 p0 c0 {3,S} {5,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-682.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,2782.5,750,1395,475,1775,1000,580,1155,1237,1373,3147,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275191,0.0883463,-0.000118429,8.11406e-08,-2.22723e-11,-81957.9,33.1246], Tmin=(100,'K'), Tmax=(886.097,'K')), NASAPolynomial(coeffs=[13.8747,0.0269566,-1.45097e-05,2.95688e-09,-2.14217e-13,-84368,-30.841], Tmin=(886.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-682.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsFHHO) + group(Cds-OdCsH) + radical(CsF1sHO2s)"""),
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
    label = 'O=CC(F)(F)OOOC(F)[C]F(6657)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {11,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
11 C u0 p0 c0 {8,D} {9,S} {14,S}
12 C u2 p0 c0 {4,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-734.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,261,493,600,1152,1365,1422,3097,2782.5,750,1395,475,1775,1000,90,150,250,247,323,377,431,433,1065,1274],'cm^-1')),
        HinderedRotor(inertia=(0.000839794,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106785,'amu*angstrom^2'), symmetry=1, barrier=(15.1883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259191,'amu*angstrom^2'), symmetry=1, barrier=(36.8189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180273,'amu*angstrom^2'), symmetry=1, barrier=(25.5931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258848,'amu*angstrom^2'), symmetry=1, barrier=(36.8306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258342,'amu*angstrom^2'), symmetry=1, barrier=(36.8024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.635643,0.108942,-0.000139233,8.84547e-08,-2.24001e-11,-88127.4,37.4801], Tmin=(100,'K'), Tmax=(958.79,'K')), NASAPolynomial(coeffs=[18.0733,0.030888,-1.71171e-05,3.54336e-09,-2.59477e-13,-91715,-51.9922], Tmin=(958.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-734.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C(F)(F)C=O(1680)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-518.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([377,639,914,1335,2782.5,750,1395,475,1775,1000,385.067,385.701,385.832,385.977],'cm^-1')),
        HinderedRotor(inertia=(0.0011382,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51134,0.0280011,-1.49649e-05,9.19102e-10,8.40317e-13,-62243.8,16.8477], Tmin=(100,'K'), Tmax=(1344.01,'K')), NASAPolynomial(coeffs=[10.7212,0.0119437,-6.39272e-06,1.30435e-09,-9.39325e-14,-65207.2,-28.0022], Tmin=(1344.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-518.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sF1sCO))"""),
)

species(
    label = 'FC1OOC1F(976)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-446.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93265,0.00400026,9.42613e-05,-1.82672e-07,1.07583e-10,-53652.7,10.67], Tmin=(10,'K'), Tmax=(549.667,'K')), NASAPolynomial(coeffs=[2.12082,0.0318962,-2.20097e-05,7.03784e-09,-8.46553e-13,-53675.7,16.3053], Tmin=(549.667,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-446.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1OOC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1O[CH]C(F)(F)OOOC1F(6658)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u0 p2 c0 {6,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
10 C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
12 C u1 p0 c0 {5,S} {11,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-926.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.414034,0.0768694,-2.76455e-05,-2.7381e-08,1.68951e-11,-111234,30.8908], Tmin=(100,'K'), Tmax=(1031.95,'K')), NASAPolynomial(coeffs=[25.2261,0.0210734,-9.90246e-06,2.08882e-09,-1.6058e-13,-118846,-104.86], Tmin=(1031.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-926.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(Cs-CsOsHH) + group(CsCFHO) + group(CsCFFO) + ring(Cyclooctane) + radical(CCsJOCs) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C1C(F)C(F)OOOC1(F)F(6659)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {7,S} {12,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u1 p2 c0 {10,S}
9  C u0 p0 c0 {1,S} {10,S} {12,S} {14,S}
10 C u0 p0 c0 {8,S} {9,S} {11,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
12 C u0 p0 c0 {4,S} {6,S} {9,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-922.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.689538,0.050106,0.000120215,-2.43025e-07,1.11001e-10,-110675,34.7101], Tmin=(100,'K'), Tmax=(914.596,'K')), NASAPolynomial(coeffs=[48.7824,-0.0274323,1.96957e-05,-3.78955e-09,2.41037e-13,-125531,-231.29], Tmin=(914.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-922.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFFO) + group(CsCFHO) + ring(Cycloheptane) + radical(CC(C)OJ) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'O=CC(F)(F)OOOC=CF(6660)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
9  C u0 p0 c0 {5,S} {11,D} {12,S}
10 C u0 p0 c0 {7,D} {8,S} {13,S}
11 C u0 p0 c0 {3,S} {9,D} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-726.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0141857,0.0925288,-0.000107331,6.31981e-08,-1.49862e-11,-87280.9,36.1877], Tmin=(100,'K'), Tmax=(1016.28,'K')), NASAPolynomial(coeffs=[15.7515,0.0305891,-1.59112e-05,3.22913e-09,-2.34404e-13,-90479.6,-39.99], Tmin=(1016.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-726.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]OOC(F)(F)C=O(6661)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u0 p2 c0 {8,D}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8 C u0 p0 c0 {5,D} {7,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-515.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,2782.5,750,1395,475,1775,1000,180,180,180,180,571.484],'cm^-1')),
        HinderedRotor(inertia=(0.00831217,'amu*angstrom^2'), symmetry=1, barrier=(1.92682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00542713,'amu*angstrom^2'), symmetry=1, barrier=(50.8779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0743114,'amu*angstrom^2'), symmetry=1, barrier=(17.2258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35371,0.0602075,-7.82013e-05,5.06303e-08,-1.2891e-11,-61864.2,25.9011], Tmin=(100,'K'), Tmax=(961.966,'K')), NASAPolynomial(coeffs=[12.3053,0.0146685,-7.19093e-06,1.41767e-09,-1.01197e-13,-63971.2,-26.5092], Tmin=(961.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = 'CHFCHF[Z](59)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-310.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([151,237,609,755,844,966,1147,1245,1323,1443,3181,3261],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383152,0.0294896,-2.94145e-05,1.64336e-08,-4.01759e-12,-36926.9,22.5083], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.34201,0.00821939,-3.17549e-06,5.49282e-10,-3.47434e-14,-38823.3,-13.1129], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-310.115,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCHF[Z]""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O=CC(F)(F)OOOC(F)=CF(6662)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {11,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {12,D}
11 C u0 p0 c0 {8,D} {9,S} {13,S}
12 C u0 p0 c0 {4,S} {10,D} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-913.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,326,540,652,719,1357,2782.5,750,1395,475,1775,1000,194,682,905,1196,1383,3221,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.572976,0.110033,-0.000162581,1.27355e-07,-4.02883e-11,-109765,39.017], Tmin=(100,'K'), Tmax=(770.605,'K')), NASAPolynomial(coeffs=[13.5287,0.0368326,-2.00894e-05,4.07839e-09,-2.93281e-13,-111938,-25.3401], Tmin=(770.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-913.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F))"""),
)

species(
    label = 'O=C[C](F)OOOC(F)[CH]F(5368)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {4,S} {10,S} {12,S}
9  C u1 p0 c0 {2,S} {5,S} {11,S}
10 C u1 p0 c0 {3,S} {8,S} {13,S}
11 C u0 p0 c0 {7,D} {9,S} {14,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-543.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,280,501,1494,1531,334,575,1197,1424,3202,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145271,0.0980656,-0.000119437,7.53171e-08,-1.92656e-11,-65279.3,36.0636], Tmin=(100,'K'), Tmax=(942.586,'K')), NASAPolynomial(coeffs=[14.8501,0.034429,-1.81651e-05,3.68876e-09,-2.67363e-13,-68106.1,-35.3936], Tmin=(942.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-543.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=CC(F)(F)OOO[CH][CH]F(6663)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {11,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
9  C u1 p0 c0 {5,S} {11,S} {12,S}
10 C u0 p0 c0 {7,D} {8,S} {13,S}
11 C u1 p0 c0 {3,S} {9,S} {14,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-512.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.059,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.826554,0.115119,-0.000175842,1.39719e-07,-4.43029e-11,-61445.5,39.1895], Tmin=(100,'K'), Tmax=(772.498,'K')), NASAPolynomial(coeffs=[14.8484,0.0339413,-1.81914e-05,3.64677e-09,-2.5968e-13,-63866.9,-32.3847], Tmin=(772.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsOsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + group(Cds-OdCsH) + radical(CCsJOO) + radical(Csj(Cs-O2sHH)(F1s)(H))"""),
)

species(
    label = '[O]OC(F)[CH]F(502)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-229.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.420111,'amu*angstrom^2'), symmetry=1, barrier=(9.65918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418454,'amu*angstrom^2'), symmetry=1, barrier=(9.62109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.65,'J/mol'), sigma=(5.51737,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.71 K, Pc=45.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66563,0.0577517,-0.000101955,9.1467e-08,-3.1227e-11,-27473.1,20.2392], Tmin=(100,'K'), Tmax=(863.103,'K')), NASAPolynomial(coeffs=[7.65514,0.0169238,-8.28529e-06,1.57108e-09,-1.06541e-13,-28020.1,-4.95473], Tmin=(863.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'O=C[C](F)F(787)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-391.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([179,346,818,1406,1524,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.55681,'amu*angstrom^2'), symmetry=1, barrier=(31.8674,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92143,0.00544306,7.55479e-05,-1.88818e-07,1.40422e-10,-47057,10.1435], Tmin=(10,'K'), Tmax=(452.346,'K')), NASAPolynomial(coeffs=[3.74847,0.0196353,-1.35046e-05,4.31321e-09,-5.19454e-13,-47170.9,9.4087], Tmin=(452.346,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-391.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OOC(F)[CH]F(4564)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {6,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-198.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,180,464.059,464.066,464.068,464.07,464.097],'cm^-1')),
        HinderedRotor(inertia=(0.0165826,'amu*angstrom^2'), symmetry=1, barrier=(2.53422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244551,'amu*angstrom^2'), symmetry=1, barrier=(37.3705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0511499,'amu*angstrom^2'), symmetry=1, barrier=(7.817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3514.65,'J/mol'), sigma=(5.71527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.98 K, Pc=42.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36452,0.0580205,-7.17947e-05,4.36937e-08,-1.03782e-11,-23756.9,24.7476], Tmin=(100,'K'), Tmax=(1033.39,'K')), NASAPolynomial(coeffs=[13.1073,0.0125657,-5.81375e-06,1.12649e-09,-8.00183e-14,-26183.8,-32.2901], Tmin=(1033.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(ROOJ) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'F[CH][CH]F(493)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-71.0739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([262,406,528,622,1148,1246,1368,1480,3164,3240,1663.85],'cm^-1')),
        HinderedRotor(inertia=(0.367109,'amu*angstrom^2'), symmetry=1, barrier=(8.44055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88778,0.0292785,-5.20594e-05,5.24809e-08,-2.00716e-11,-8512.84,14.1703], Tmin=(100,'K'), Tmax=(832.93,'K')), NASAPolynomial(coeffs=[3.49548,0.0155532,-7.88013e-06,1.54339e-09,-1.07576e-13,-8239.2,13.6], Tmin=(832.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.0739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sH)"""),
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
    label = 'F[CH]C(F)OOO[C](F)F(5935)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
9  C u1 p0 c0 {2,S} {8,S} {12,S}
10 C u1 p0 c0 {3,S} {4,S} {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-579.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,493,600,700,1144,1293,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0391999,0.0898682,-0.000112103,6.84585e-08,-1.64544e-11,-69516.6,34.0646], Tmin=(100,'K'), Tmax=(1015.43,'K')), NASAPolynomial(coeffs=[17.3555,0.0216562,-1.13414e-05,2.30532e-09,-1.6759e-13,-73033.4,-49.7421], Tmin=(1015.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFHO) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'O=CC(F)(F)OOO[C](F)[CH]F(6664)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {11,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u1 p0 c0 {3,S} {6,S} {12,S}
11 C u0 p0 c0 {8,D} {9,S} {13,S}
12 C u1 p0 c0 {4,S} {10,S} {14,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {12,S}
"""),
    E0 = (-711.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,395,473,707,1436,2782.5,750,1395,475,1775,1000,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34602,0.127419,-0.000202964,1.63962e-07,-5.23763e-11,-85384.6,41.5545], Tmin=(100,'K'), Tmax=(769.124,'K')), NASAPolynomial(coeffs=[16.8883,0.0325826,-1.7999e-05,3.62983e-09,-2.58566e-13,-88189.3,-41.628], Tmin=(769.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'O=[C]C(F)(F)OOOC(F)[CH]F(6665)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {12,D}
9  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {6,S} {12,S}
11 C u1 p0 c0 {4,S} {9,S} {14,S}
12 C u1 p0 c0 {8,D} {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-735.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,254,363,525,544,658,1137,1402,334,575,1197,1424,3202,1855,455,950,180,180,180,330.332,1426.39,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.144962,'amu*angstrom^2'), symmetry=1, barrier=(3.33296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144962,'amu*angstrom^2'), symmetry=1, barrier=(3.33296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144962,'amu*angstrom^2'), symmetry=1, barrier=(3.33296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144962,'amu*angstrom^2'), symmetry=1, barrier=(3.33296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144962,'amu*angstrom^2'), symmetry=1, barrier=(3.33296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144962,'amu*angstrom^2'), symmetry=1, barrier=(3.33296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95832,0.149077,-0.000282492,2.56117e-07,-8.71683e-11,-88215.2,42.247], Tmin=(100,'K'), Tmax=(872.526,'K')), NASAPolynomial(coeffs=[14.678,0.0341327,-1.83932e-05,3.52373e-09,-2.37602e-13,-89646,-27.3085], Tmin=(872.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-735.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-F1sO2sH)(F1s)(H)) + radical(COj(Cs-F1sF1sO2s)(O2d))"""),
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
    label = 'O=C=C(F)OOOC(F)[CH]F(6666)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {6,S} {10,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {11,D}
8  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
9  C u1 p0 c0 {2,S} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {5,S} {11,D}
11 C u0 p0 c0 {7,D} {10,D}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
"""),
    E0 = (-460.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,334,575,1197,1424,3202,197,221,431,657,2120,512.5,787.5,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.267577,0.102029,-0.000145034,1.06454e-07,-3.13309e-11,-55244.9,37.5414], Tmin=(100,'K'), Tmax=(828.259,'K')), NASAPolynomial(coeffs=[14.3612,0.0313833,-1.70985e-05,3.48193e-09,-2.51405e-13,-57668.3,-30.2782], Tmin=(828.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-460.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + missing(O2d-Cdd) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
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
    label = 'O=CC(F)(F)OOOC(F)[C]F-2(6667)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {11,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
11 C u0 p0 c0 {8,D} {9,S} {14,S}
12 C u0 p1 c0 {4,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
"""),
    E0 = (-752.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,617,898,1187,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (190.05,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.5771,0.107498,-0.000138147,8.88009e-08,-2.2765e-11,-90385.6,38.593], Tmin=(100,'K'), Tmax=(947.657,'K')), NASAPolynomial(coeffs=[17.581,0.0308563,-1.68369e-05,3.4629e-09,-2.52659e-13,-93827.2,-48.0336], Tmin=(947.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-752.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFHO) + group(Cds-OdCsH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'O=CC(F)(F)OOO[C](F)CF(6668)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {11,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {12,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {3,S} {11,S} {13,S} {14,S}
11 C u1 p0 c0 {4,S} {6,S} {10,S}
12 C u0 p0 c0 {8,D} {9,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-910.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.882078,0.116281,-0.000159892,1.13223e-07,-3.2161e-11,-109291,39.554], Tmin=(100,'K'), Tmax=(857.046,'K')), NASAPolynomial(coeffs=[16.314,0.036021,-1.94182e-05,3.95054e-09,-2.85573e-13,-112238,-40.7541], Tmin=(857.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-910.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)(F)OOOC(F)CF(6669)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {11,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {12,D}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
12 C u1 p0 c0 {8,D} {11,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
"""),
    E0 = (-933.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7514,0.141226,-0.000252236,2.24197e-07,-7.613e-11,-112111,41.1525], Tmin=(100,'K'), Tmax=(853.339,'K')), NASAPolynomial(coeffs=[14.5439,0.0367903,-1.93488e-05,3.7326e-09,-2.55188e-13,-113870,-28.8941], Tmin=(853.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-933.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(COj(Cs-F1sF1sO2s)(O2d))"""),
)

species(
    label = 'O=C[C](F)OOOC(F)C(F)F(6670)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {11,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {12,D}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
10 C u0 p0 c0 {2,S} {3,S} {9,S} {14,S}
11 C u1 p0 c0 {4,S} {6,S} {12,S}
12 C u0 p0 c0 {8,D} {11,S} {15,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-958.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.325347,0.101018,-0.000113346,6.48511e-08,-1.50604e-11,-115176,37.2375], Tmin=(100,'K'), Tmax=(1031.71,'K')), NASAPolynomial(coeffs=[16.4218,0.0360888,-1.89468e-05,3.85369e-09,-2.79946e-13,-118632,-44.0808], Tmin=(1031.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-958.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=CC(F)(F)OOO[CH]C(F)F(5996)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {7,S} {11,S}
7  O u0 p2 c0 {5,S} {6,S}
8  O u0 p2 c0 {12,D}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {11,S} {13,S}
11 C u1 p0 c0 {6,S} {10,S} {14,S}
12 C u0 p0 c0 {8,D} {9,S} {15,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {12,S}
"""),
    E0 = (-927.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (191.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.929582,0.117484,-0.000164698,1.19694e-07,-3.49028e-11,-111403,39.39], Tmin=(100,'K'), Tmax=(835.716,'K')), NASAPolynomial(coeffs=[15.9569,0.0366612,-1.96325e-05,3.9741e-09,-2.86137e-13,-114226,-39.0471], Tmin=(835.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-927.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsOsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CCsJOO)"""),
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
    E0 = (-441.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-202.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (192.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (76.0447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-68.6614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-8.93432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-67.1352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-121.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-415.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-454.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-435.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-208.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-414.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-301.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-70.6217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-38.9285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-346.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-189.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-384.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-185.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-146.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-99.2044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-119.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-151.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-143.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-132.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-361.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-420.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-333.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-327.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['[O]OC(F)(F)C=O(1651)', 'FC1OC1F(977)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.38e+12,'s^-1','*|/',1.2), n=0, Ea=(63.953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_S;C_sec_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'F[CH]C(F)OOOC(F)F(5928)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.81927e-16,'m^3/(mol*s)'), n=6.80628, Ea=(313.044,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C[C]F-2(1630)', 'F[CH]C(F)OOOF(5927)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(11.5935,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[CH]C(F)OOOC(F)=COF(6654)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(130.027,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(F)OOOOC=C(F)F(6655)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(97.3416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[CH]C(F)OOOC(F)(F)C=O(6656)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]F(223)', 'O=CC(F)(F)OOO[CH]F(3089)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', 'O=CC(F)(F)OOOC(F)[C]F(6657)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Hrad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['[O]C(F)(F)C=O(1680)', 'FC1OOC1F(976)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0, Ea=(89.9403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO;C_sec_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['FC1O[CH]C(F)(F)OOOC1F(6658)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(50.8858,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['[O]C1C(F)C(F)OOOC1(F)F(6659)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19e+11,'s^-1'), n=0.08, Ea=(69.8728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8;multiplebond_intra;radadd_intra_cs] for rate rule [R8;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'O=CC(F)(F)OOOC=CF(6660)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(45.2717,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OOC(F)(F)C=O(6661)', 'CHFCHF[Z](59)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000422724,'m^3/(mol*s)'), n=2.88896, Ea=(10.4335,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.38434251561780236, var=0.3755598567717395, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0_N-4R!H->S_3OS->O_Ext-2CNS-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_4R!H-u0_N-4R!H->S_3OS->O_Ext-2CNS-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'O=CC(F)(F)OOOC(F)=CF(6662)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(21820,'m^3/(mol*s)'), n=0.859, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_4R!H->O"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'O=C[C](F)OOOC(F)[CH]F(5368)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -26.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=CC(F)(F)OOO[CH][CH]F(6663)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -34.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)(F)C=O(1680)', '[O]OC(F)[CH]F(502)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=C[C](F)F(787)', '[O]OOC(F)[CH]F(4564)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OC(F)(F)C=O(1651)', '[O]C(F)[CH]F(975)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(2.77885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OOC(F)(F)C=O(6661)', 'F[CH][CH]F(493)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HCO(15)', 'F[CH]C(F)OOO[C](F)F(5935)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -5.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', 'O=CC(F)(F)OOO[C](F)[CH]F(6664)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R
Ea raised from -0.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'O=[C]C(F)(F)OOOC(F)[CH]F(6665)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.7313e+28,'m^3/(mol*s)'), n=-7.39166, Ea=(3.45096,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3166437201132897, var=29.409938925631906, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=C=C(F)OOOC(F)[CH]F(6666)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(190.048,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CHF(40)', 'O=CC(F)(F)OOO[CH]F(3089)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(5)', 'O=CC(F)(F)OOOC(F)[C]F-2(6667)'],
    products = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.75,'m^3/(mol*s)'), n=-0.32, Ea=(8.03896,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_3R->H_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['O=CC(F)(F)OOO[C](F)CF(6668)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(34257.4,'s^-1'), n=2.45724, Ea=(144.304,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out] for rate rule [R2H_S;C_rad_out_1H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['O=[C]C(F)(F)OOOC(F)CF(6669)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(383890,'s^-1'), n=1.48889, Ea=(84.9143,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_single;XH_out] for rate rule [R7H;C_rad_out_1H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['O=C[C](F)OOOC(F)C(F)F(6670)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.29706e+07,'s^-1'), n=1.15307, Ea=(171.976,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=CC(F)(F)OOOC(F)[CH]F(5994)'],
    products = ['O=CC(F)(F)OOO[CH]C(F)F(5996)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(177.629,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #2266',
    isomers = [
        'O=CC(F)(F)OOOC(F)[CH]F(5994)',
    ],
    reactants = [
        ('[O]OC(F)(F)C=O(1651)', 'FC1OC1F(977)'),
        ('[O]OC(F)(F)C=O(1651)', '[O]C(F)[CH]F(975)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2266',
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

