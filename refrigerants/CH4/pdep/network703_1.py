species(
    label = 'O=[C]OO[C]=C=O(1309)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {7,D}
6 C u1 p0 c0 {2,S} {3,D}
7 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (136.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,1855,455,950,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(1.80241,'amu*angstrom^2'), symmetry=1, barrier=(41.4409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80302,'amu*angstrom^2'), symmetry=1, barrier=(41.4549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80357,'amu*angstrom^2'), symmetry=1, barrier=(41.4675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6365,0.0484287,-5.56198e-05,2.99527e-08,-6.20171e-12,16545.9,24.96], Tmin=(100,'K'), Tmax=(1189.84,'K')), NASAPolynomial(coeffs=[13.9689,0.00696951,-3.35292e-06,6.67344e-10,-4.84501e-14,13611.2,-36.6799], Tmin=(1189.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(Cds-OdOsH) + missing(Cdd-CdO2d) + radical(C=CJO) + radical((O)CJOC)"""),
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
    label = 'O=C=C=O(813)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (63.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2110,2130,495,530,650,925,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2387.04,'J/mol'), sigma=(4.99307,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=372.85 K, Pc=43.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5355,0.024033,-4.26139e-05,3.87325e-08,-1.34306e-11,7697.1,25.5769], Tmin=(100,'K'), Tmax=(857.394,'K')), NASAPolynomial(coeffs=[4.78324,0.00773532,-3.93446e-06,7.52107e-10,-5.12483e-14,7525.26,16.3243], Tmin=(857.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.731,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""OCCO(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O][C]=O(213)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3422.77,'J/mol'), sigma=(5.42046,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=534.63 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186753,1.03202e-05,-1.52373e-08,5.80537e-12,3804.5,8.40407], Tmin=(100,'K'), Tmax=(1021.26,'K')), NASAPolynomial(coeffs=[6.36178,0.000422737,-4.06602e-07,1.525e-10,-1.51974e-14,2816.75,-6.43921], Tmin=(1021.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[O]O[C]=O(131)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {1,S} {3,D}
"""),
    E0 = (106.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,262.32,262.473],'cm^-1')),
        HinderedRotor(inertia=(0.00244921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0088,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3570.08,'J/mol'), sigma=(5.61676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.64 K, Pc=45.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24172,0.0164884,-1.87037e-05,1.04488e-08,-2.29625e-12,12889.7,11.1068], Tmin=(100,'K'), Tmax=(1108.87,'K')), NASAPolynomial(coeffs=[6.67162,0.0041158,-1.96699e-06,3.86455e-10,-2.76579e-14,12129,-5.79494], Tmin=(1108.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""CO3t2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[C]=C=O(1258)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 C u2 p0 c0 {2,D}
"""),
    E0 = (368.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,444.973],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86278,0.0119701,-1.80851e-05,1.52778e-08,-5.20063e-12,44312.6,8.89759], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.42468,0.00185394,-5.17933e-07,6.77646e-11,-3.53315e-15,43716.1,-3.69608], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(368.265,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[C]=O(117)',
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
    label = '[O]O[C]=C=O(1312)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
"""),
    E0 = (322.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1685,370,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.243089,'amu*angstrom^2'), symmetry=1, barrier=(5.58909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72851,0.0348827,-7.37524e-05,7.20645e-08,-2.53326e-11,38801.1,18.3764], Tmin=(100,'K'), Tmax=(910.539,'K')), NASAPolynomial(coeffs=[4.31527,0.0105752,-5.14874e-06,9.24356e-10,-5.88728e-14,39230.8,14.8163], Tmin=(910.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'O=C=C1OOC1=O(1310)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {6,S} {7,D}
6 C u0 p0 c0 {2,S} {3,D} {5,S}
7 C u0 p0 c0 {4,D} {5,D}
"""),
    E0 = (-104.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35714,0.0233919,1.63797e-05,-4.48377e-08,2.01082e-11,-12478.4,20.6892], Tmin=(100,'K'), Tmax=(971.382,'K')), NASAPolynomial(coeffs=[14.7526,0.00357035,-1.2233e-06,3.31205e-10,-3.2363e-14,-16359.5,-46.3338], Tmin=(971.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-O2d(Cds-Cds)O2s) + missing(Cdd-CdO2d) + ring(Cyclobutane)"""),
)

species(
    label = 'O=[C]OOC1=[C]O1(1313)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u1 p0 c0 {1,S} {5,D}
7 C u1 p0 c0 {3,S} {4,D}
"""),
    E0 = (433.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51929,0.0647821,-0.000127658,1.25002e-07,-4.60034e-11,52196,24.9091], Tmin=(100,'K'), Tmax=(839.657,'K')), NASAPolynomial(coeffs=[5.62212,0.0229999,-1.32912e-05,2.6565e-09,-1.8538e-13,52290.9,10.5002], Tmin=(839.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdOsH) + ring(oxirene) + radical(C=CJO) + radical((O)CJOC)"""),
)

species(
    label = 'O=C1O[C]=[C]OO1(1314)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 C u1 p0 c0 {2,S} {7,D}
7 C u1 p0 c0 {3,S} {6,D}
"""),
    E0 = (22.2558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52712,0.00695283,9.0001e-05,-1.30019e-07,5.106e-11,2752.74,19.3898], Tmin=(100,'K'), Tmax=(971.462,'K')), NASAPolynomial(coeffs=[19.5034,0.00138695,-7.41211e-07,4.22147e-10,-5.14722e-14,-3581.36,-77.6442], Tmin=(971.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.2558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + ring(Cyclohexane) + radical(C=CJO) + radical(C=CJO)"""),
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
    label = '[O]C#C[O](814)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
"""),
    E0 = (9.55717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,344.097,344.099,344.1],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09811,0.0228497,-4.17927e-05,3.86217e-08,-1.34608e-11,1179.06,7.869], Tmin=(100,'K'), Tmax=(871.084,'K')), NASAPolynomial(coeffs=[4.9831,0.00749284,-3.80902e-06,7.20359e-10,-4.85535e-14,1104.9,0.494495], Tmin=(871.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.55717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""OCCO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O(7)',
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
    label = '[C]#COO[C]=O(1315)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {2,S} {3,D}
5 C u0 p0 c0 {1,S} {6,T}
6 C u1 p0 c0 {5,T}
"""),
    E0 = (493.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1855,455,950,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.90054,'amu*angstrom^2'), symmetry=1, barrier=(43.6972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90026,'amu*angstrom^2'), symmetry=1, barrier=(43.6908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9029,'amu*angstrom^2'), symmetry=1, barrier=(43.7514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0301,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55797,0.0466397,-5.49943e-05,2.96166e-08,-5.99943e-12,59396.7,19.8005], Tmin=(100,'K'), Tmax=(1285.58,'K')), NASAPolynomial(coeffs=[15.2327,0.00198959,-4.44474e-07,5.67263e-11,-3.7399e-15,56054.4,-48.9319], Tmin=(1285.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(120.56,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsCt) + group(Cds-OdOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical((O)CJOC) + radical(Acetyl)"""),
)

species(
    label = 'O=C1[C]OOC1=O(1316)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 C u2 p0 c0 {2,S} {6,S}
"""),
    E0 = (27.3355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26892,0.0146154,6.78567e-05,-1.11118e-07,4.60375e-11,3371.4,17.8242], Tmin=(100,'K'), Tmax=(954.955,'K')), NASAPolynomial(coeffs=[20.5426,-0.00269597,2.01116e-06,-1.99192e-10,-3.69736e-15,-2819.45,-83.6341], Tmin=(954.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.3355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + ring(Cyclopentane) + radical(CH2_triplet)"""),
)

species(
    label = 'O=[C][C]1OOC1=O(1317)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {3,D} {5,S}
7 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-35.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40742,0.0219107,2.38581e-05,-5.0127e-08,2.09311e-11,-4226.63,20.6568], Tmin=(100,'K'), Tmax=(992.551,'K')), NASAPolynomial(coeffs=[13.7037,0.00830261,-3.81024e-06,8.54055e-10,-6.98359e-14,-8041.19,-41.6765], Tmin=(992.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(C2CsJOO) + radical(C=OCCJ=O)"""),
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
    E0 = (22.5682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (360.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (647.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (30.8525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (320.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (88.8342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (100.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (60.2076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (22.5682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (40.9172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (621.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (81.1442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (103.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]OO[C]=C=O(1309)'],
    products = ['CO2(15)', 'O=C=C=O(813)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]O[C]=O(131)', '[C]=C=O(1258)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(117)', '[O]O[C]=C=O(1312)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]OO[C]=C=O(1309)'],
    products = ['O=C=C1OOC1=O(1310)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C]OO[C]=C=O(1309)'],
    products = ['O=[C]OOC1=[C]O1(1313)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.631e+12,'s^-1'), n=0.216, Ea=(298.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3;multiplebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]OO[C]=C=O(1309)'],
    products = ['O=C1O[C]=[C]OO1(1314)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.32511e+10,'s^-1'), n=0.274726, Ea=(66.2661,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;multiplebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CO(14)', '[O]O[C]=C=O(1312)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(34.1,'m^3/(mol*s)'), n=8.73864e-09, Ea=(11.4672,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]=O(213)', 'O=C=C=O(813)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.64086e-07,'m^3/(mol*s)'), n=3.71185, Ea=(79.2107,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CO2(15)', '[O]C#C[O](814)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.64086e-07,'m^3/(mol*s)'), n=3.71185, Ea=(530.418,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 4.0
Ea raised from 526.0 to 530.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=O(213)', '[O]C#C[O](814)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.95326e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(114.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O(7)', '[C]#COO[C]=O(1315)'],
    products = ['O=[C]OO[C]=C=O(1309)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Ct_rad;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C]OO[C]=C=O(1309)'],
    products = ['O=C1[C]OOC1=O(1316)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra;radadd_intra] for rate rule [R5_SS_T;triplebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]OO[C]=C=O(1309)'],
    products = ['O=[C][C]1OOC1=O(1317)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_T;triplebond_intra;radadd_intra_CO]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = 'PDepNetwork #703',
    isomers = [
        'O=[C]OO[C]=C=O(1309)',
    ],
    reactants = [
        ('CO2(15)', 'O=C=C=O(813)'),
        ('[O][C]=O(213)', 'O=C=C=O(813)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #703',
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

