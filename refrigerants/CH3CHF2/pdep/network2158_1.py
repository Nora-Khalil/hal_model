species(
    label = 'F[C](F)OC1O[C]1F(5674)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {6,S} {8,S}
6 C u0 p0 c0 {4,S} {5,S} {7,S} {9,S}
7 C u1 p0 c0 {1,S} {4,S} {6,S}
8 C u1 p0 c0 {2,S} {3,S} {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-460.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,493,600,700,1144,1293,180,180,180,942.709,942.722,942.831,942.833,943.013],'cm^-1')),
        HinderedRotor(inertia=(0.133603,'amu*angstrom^2'), symmetry=1, barrier=(3.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133748,'amu*angstrom^2'), symmetry=1, barrier=(3.07512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75354,0.0545273,-6.76791e-05,4.58918e-08,-1.29422e-11,-55349.1,23.0283], Tmin=(100,'K'), Tmax=(849.595,'K')), NASAPolynomial(coeffs=[8.4013,0.0232293,-1.24221e-05,2.53307e-09,-1.83774e-13,-56478.7,-7.9599], Tmin=(849.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-460.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + ring(Cs(O2)-O2s-Cs(F)) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'CF2O(49)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-618.61,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=CO1(4047)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (74.4547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,326,540,652,719,1357,615.392,2893.28],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3010.76,'J/mol'), sigma=(4.86933,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=470.27 K, Pc=59.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.102,0.0258723,-5.29144e-05,5.76368e-08,-2.28587e-11,8981.24,10.2596], Tmin=(100,'K'), Tmax=(846.704,'K')), NASAPolynomial(coeffs=[1.99159,0.016003,-8.6525e-06,1.7026e-09,-1.18158e-13,9711.08,18.6315], Tmin=(846.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.4547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + ring(oxirene)"""),
)

species(
    label = 'F[C]F(2927)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = '[O]C1O[C]1F(5677)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-20.9123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,395,473,707,1436,180,1138.2,1138.29,1138.41,1138.47,1138.68],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22088,0.0176379,-7.66961e-06,4.58929e-10,2.26595e-13,-2487.61,13.5555], Tmin=(100,'K'), Tmax=(1758.17,'K')), NASAPolynomial(coeffs=[9.09322,0.00918078,-4.63744e-06,8.95342e-10,-6.10015e-14,-5310.32,-20.2438], Tmin=(1758.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.9123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'FC1(F)OC2OC21F(5676)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {7,S} {8,S}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-782.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51195,0.045121,-2.71472e-05,4.99751e-10,2.77415e-12,-94071.9,16.2913], Tmin=(100,'K'), Tmax=(1165.73,'K')), NASAPolynomial(coeffs=[15.8316,0.0128522,-7.32849e-06,1.57741e-09,-1.18755e-13,-98556.5,-59.9041], Tmin=(1165.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-782.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsOsH) + group(CsCFFO) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC1=C(OC(F)F)O1(5703)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {7,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
7 C u0 p0 c0 {4,S} {5,S} {8,D}
8 C u0 p0 c0 {3,S} {5,S} {7,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-472.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968921,0.0866084,-0.00019148,2.0112e-07,-7.63681e-11,-56778.4,22.1255], Tmin=(100,'K'), Tmax=(861.293,'K')), NASAPolynomial(coeffs=[1.18205,0.0372577,-2.13084e-05,4.21033e-09,-2.90151e-13,-55021.3,31.5423], Tmin=(861.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsFFHO) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(oxirene)"""),
)

species(
    label = 'O=C(F)[CH]O[C](F)F(5704)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {4,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 C u1 p0 c0 {2,S} {3,S} {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-610.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,493,600,700,1144,1293,244.153,244.166,244.198,244.198],'cm^-1')),
        HinderedRotor(inertia=(0.237006,'amu*angstrom^2'), symmetry=1, barrier=(10.0272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00235805,'amu*angstrom^2'), symmetry=1, barrier=(10.0269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645282,'amu*angstrom^2'), symmetry=1, barrier=(27.298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888445,0.0741621,-0.000124068,1.048e-07,-3.46432e-11,-73370.2,26.6211], Tmin=(100,'K'), Tmax=(800.21,'K')), NASAPolynomial(coeffs=[11.1627,0.0179605,-9.63843e-06,1.90235e-09,-1.32985e-13,-74859.4,-19.6874], Tmin=(800.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-610.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFHO) + group(COCsFO) + radical(CCsJOCs) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'F[C]1[CH]O1(5679)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18921,0.0134531,-5.29457e-07,-1.06374e-08,5.34941e-12,11790.3,13.5601], Tmin=(100,'K'), Tmax=(993.343,'K')), NASAPolynomial(coeffs=[8.19763,0.00359113,-1.19992e-06,2.57126e-10,-2.11232e-14,10286.8,-13.1285], Tmin=(993.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.7552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O][C](F)F(4778)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-255.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53788e-05,4.83907e-09,3.86802e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58999,0.00108969,-4.53602e-07,1.24872e-10,-1.13713e-14,-32130.4,-16.3888], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
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
    label = 'F[C](F)OC1=C(F)O1(5705)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {6,S} {8,S}
6 C u0 p0 c0 {4,S} {5,S} {7,D}
7 C u0 p0 c0 {1,S} {4,S} {6,D}
8 C u1 p0 c0 {2,S} {3,S} {5,S}
"""),
    E0 = (-258.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,493,600,700,1144,1293,180,180,180,180,2129.7,2132.9],'cm^-1')),
        HinderedRotor(inertia=(0.161435,'amu*angstrom^2'), symmetry=1, barrier=(3.7117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161102,'amu*angstrom^2'), symmetry=1, barrier=(3.70406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21593,0.0800985,-0.000180183,1.91056e-07,-7.3184e-11,-30999.7,23.5129], Tmin=(100,'K'), Tmax=(855.484,'K')), NASAPolynomial(coeffs=[1.13831,0.0347285,-2.04436e-05,4.08405e-09,-2.83378e-13,-29312.9,33.6562], Tmin=(855.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsFFHO) + group(Cds-CdsCsCs) + group(CdCFO) + ring(oxirene) + radical(Csj(F1s)(F1s)(O2s-Cd)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    label = 'F[C](F)O[C]1OC1F(5706)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {7,S} {8,S}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7 C u1 p0 c0 {4,S} {5,S} {6,S}
8 C u1 p0 c0 {2,S} {3,S} {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-480.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94452,0.0507171,-5.92532e-05,3.90834e-08,-1.10338e-11,-57756.8,23.2175], Tmin=(100,'K'), Tmax=(837.353,'K')), NASAPolynomial(coeffs=[7.26849,0.025285,-1.36958e-05,2.81283e-09,-2.05043e-13,-58648.5,-1.52262], Tmin=(837.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-480.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + ring(Cs(O2)-O2s-Cs(F)) + radical(Cs_P) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'F[C]1O[C]1OC(F)F(5707)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {7,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
7 C u1 p0 c0 {4,S} {5,S} {8,S}
8 C u1 p0 c0 {3,S} {5,S} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-467.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77395,0.0530777,-6.02963e-05,3.5832e-08,-8.75579e-12,-56138.2,23.1719], Tmin=(100,'K'), Tmax=(976.794,'K')), NASAPolynomial(coeffs=[9.59477,0.021051,-1.11146e-05,2.26504e-09,-1.64628e-13,-57666,-14.3753], Tmin=(976.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-467.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + ring(Cs(O2)-O2s-Cs(F)) + radical(Cs_P) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring)"""),
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
    E0 = (-224.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (248.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-216.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-161.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-188.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-162.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (55.0578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (189.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (78.3583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (18.5926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-50.5771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-83.7288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)OC1O[C]1F(5674)'],
    products = ['CF2O(49)', 'FC1=CO1(4047)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C]F(2927)', '[O]C1O[C]1F(5677)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C](F)OC1O[C]1F(5674)'],
    products = ['FC1(F)OC2OC21F(5676)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)OC1O[C]1F(5674)'],
    products = ['FC1=C(OC(F)F)O1(5703)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)[CH]O[C](F)F(5704)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CF2O(49)', 'F[C]1[CH]O1(5679)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.66944e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(122.245,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C](F)F(4778)', 'FC1=CO1(4047)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', 'F[C](F)OC1=C(F)O1(5705)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C](F)F(4778)', 'F[C]1[CH]O1(5679)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CF2(43)', '[O]C1O[C]1F(5677)'],
    products = ['F[C](F)OC1O[C]1F(5674)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(7.13682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[C](F)OC1O[C]1F(5674)'],
    products = ['F[C](F)O[C]1OC1F(5706)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.06381e+10,'s^-1'), n=0.647667, Ea=(174.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_NonDe] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C](F)OC1O[C]1F(5674)'],
    products = ['F[C]1O[C]1OC(F)F(5707)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.95012e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2158',
    isomers = [
        'F[C](F)OC1O[C]1F(5674)',
    ],
    reactants = [
        ('CF2O(49)', 'FC1=CO1(4047)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2158',
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

