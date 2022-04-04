species(
    label = 'O=COOC(=O)F(8762)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {6,D}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {4,D} {8,S}
7 C u0 p0 c0 {1,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-829.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(1.83751,'amu*angstrom^2'), symmetry=1, barrier=(42.2479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83581,'amu*angstrom^2'), symmetry=1, barrier=(42.2088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83805,'amu*angstrom^2'), symmetry=1, barrier=(42.2604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8989,0.0221911,5.87633e-05,-1.06919e-07,4.56845e-11,-99707.2,22.3111], Tmin=(100,'K'), Tmax=(950.526,'K')), NASAPolynomial(coeffs=[22.1967,-0.00264097,2.34275e-06,-2.92015e-10,3.71291e-15,-106303,-88.9808], Tmin=(950.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-829.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(Cds-OdOsH) + group(COFOO)"""),
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
    label = 'O=C(F)OO(2257)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-521.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(0.0284926,'amu*angstrom^2'), symmetry=1, barrier=(17.028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.61347,'amu*angstrom^2'), symmetry=1, barrier=(60.0888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92212,0.00493879,7.38492e-05,-1.69733e-07,1.14297e-10,-62737,9.3347], Tmin=(10,'K'), Tmax=(510.617,'K')), NASAPolynomial(coeffs=[4.17798,0.0187211,-1.30128e-05,4.22274e-09,-5.16224e-13,-62969,6.25709], Tmin=(510.617,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-521.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC(F)OO""", comment="""Thermo library: 2-BTP_G4"""),
)

species(
    label = 'O=COOF(8964)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {2,S} {4,D} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-282.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,2782.5,750,1395,475,1775,1000,405.417],'cm^-1')),
        HinderedRotor(inertia=(0.388749,'amu*angstrom^2'), symmetry=1, barrier=(44.8739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383992,'amu*angstrom^2'), symmetry=1, barrier=(44.9232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73077,0.015503,2.68901e-05,-5.20542e-08,2.2149e-11,-33863.6,15.2257], Tmin=(100,'K'), Tmax=(964.099,'K')), NASAPolynomial(coeffs=[13.3654,0.00228028,-6.12678e-07,2.07508e-10,-2.33713e-14,-37350.2,-43.1389], Tmin=(964.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2sFO) + group(Cds-OdOsH)"""),
)

species(
    label = 'F[C]1OO[CH]OO1(8965)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u0 p2 c0 {4,S} {7,S}
6 C u1 p0 c0 {2,S} {4,S} {8,S}
7 C u1 p0 c0 {1,S} {3,S} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-80.9301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,482,586,761,1411,435.418,435.419,435.42,435.42,435.42,435.42,435.42,435.421,435.421,435.421,435.421,435.422],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85668,0.0344762,5.37336e-06,-3.11368e-08,1.29139e-11,-9645.78,15.5131], Tmin=(100,'K'), Tmax=(1098.31,'K')), NASAPolynomial(coeffs=[14.9372,0.0160655,-9.39992e-06,2.05983e-09,-1.57494e-13,-14281.9,-56.8448], Tmin=(1098.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.9301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-OsOsHH) + group(CsFHOO) + ring(Cyclohexane) + radical(OCJO) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]C(F)OO[C]=O(8800)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7 C u1 p0 c0 {3,S} {5,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-353.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,451,553,637,1069,1180,1265,1301,3056,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.85382,'amu*angstrom^2'), symmetry=1, barrier=(42.623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85275,'amu*angstrom^2'), symmetry=1, barrier=(42.5983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85249,'amu*angstrom^2'), symmetry=1, barrier=(42.5925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49248,0.0444225,-2.41231e-05,-9.35142e-09,8.59782e-12,-42382.8,23.1362], Tmin=(100,'K'), Tmax=(991.809,'K')), NASAPolynomial(coeffs=[16.9024,0.006245,-2.63786e-06,5.75716e-10,-4.70108e-14,-46618.6,-57.0249], Tmin=(991.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(Cds-OdOsH) + radical(O2sj(Cs-F1sO2sH)) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OO[C](O)F(8760)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {2,S} {4,S}
7 C u1 p0 c0 {3,S} {5,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-402.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,482,586,761,1411,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.82132,'amu*angstrom^2'), symmetry=1, barrier=(41.8758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82196,'amu*angstrom^2'), symmetry=1, barrier=(41.8904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82175,'amu*angstrom^2'), symmetry=1, barrier=(41.8857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82218,'amu*angstrom^2'), symmetry=1, barrier=(41.8955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3872.75,'J/mol'), sigma=(6.04376,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.91 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30622,0.0504256,-4.2566e-05,1.01026e-08,1.64728e-12,-48312.7,24.4522], Tmin=(100,'K'), Tmax=(1020.88,'K')), NASAPolynomial(coeffs=[16.9186,0.00635124,-2.92872e-06,6.23603e-10,-4.88818e-14,-52391.3,-55.5548], Tmin=(1020.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(Cds-OdOsH) + radical(CsF1sO2sO2s) + radical((O)CJOC)"""),
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
    label = 'O=[C]OOC=O(8966)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,D} {7,S}
6 C u1 p0 c0 {2,S} {4,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-410.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.89171,'amu*angstrom^2'), symmetry=1, barrier=(43.4942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89226,'amu*angstrom^2'), symmetry=1, barrier=(43.5068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0268,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26635,0.0157923,6.27644e-05,-1.06118e-07,4.47151e-11,-49321.7,18.9319], Tmin=(100,'K'), Tmax=(945.571,'K')), NASAPolynomial(coeffs=[20.082,-0.00293736,2.63337e-06,-3.79914e-10,1.17882e-14,-55222.8,-79.4094], Tmin=(945.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(Cds-OdOsH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C=O(322)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-138.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51545,0.0112656,-9.58083e-06,4.36318e-09,-8.44581e-13,-16672.2,7.37034], Tmin=(100,'K'), Tmax=(1184.07,'K')), NASAPolynomial(coeffs=[5.09941,0.00591472,-2.80224e-06,5.4664e-10,-3.87737e-14,-17047.4,-0.538979], Tmin=(1184.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""formyloxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]OC(=O)F(277)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-339.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(1.35484,'amu*angstrom^2'), symmetry=1, barrier=(31.1504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92964,0.00448623,6.20575e-05,-1.45434e-07,9.86982e-11,-40867.8,10.2015], Tmin=(10,'K'), Tmax=(509.958,'K')), NASAPolynomial(coeffs=[4.25304,0.0157357,-1.15832e-05,3.84885e-09,-4.74388e-13,-41080,7.10146], Tmin=(509.958,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-339.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]OC(DO)F""", comment="""Thermo library: 2-BTP_G4"""),
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
    label = '[O]OC=O(320)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-118.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,524.127,524.138],'cm^-1')),
        HinderedRotor(inertia=(0.141181,'amu*angstrom^2'), symmetry=1, barrier=(27.5177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3570.08,'J/mol'), sigma=(5.61676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.64 K, Pc=45.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.11871,0.0145187,-1.84988e-07,-1.22536e-08,6.08153e-12,-14223.7,11.3177], Tmin=(100,'K'), Tmax=(990.678,'K')), NASAPolynomial(coeffs=[8.52145,0.00409734,-1.65621e-06,3.44828e-10,-2.71434e-14,-15853.2,-17.5185], Tmin=(990.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formylperoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'O=[C]OOC(=O)F(8799)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {6,D}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 C u1 p0 c0 {3,S} {5,D}
"""),
    E0 = (-633.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,482,664,788,1296,1923,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.84063,'amu*angstrom^2'), symmetry=1, barrier=(42.3198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84232,'amu*angstrom^2'), symmetry=1, barrier=(42.3585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84066,'amu*angstrom^2'), symmetry=1, barrier=(42.3203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.017,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69273,0.029917,3.07951e-05,-7.8832e-08,3.6421e-11,-76075.3,21.6785], Tmin=(100,'K'), Tmax=(943.488,'K')), NASAPolynomial(coeffs=[22.4991,-0.00567578,3.72802e-06,-5.96522e-10,2.80335e-14,-82343.4,-89.9009], Tmin=(943.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-633.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-O2s(Cds-O2d)) + group(COFOO) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    E0 = (-497.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-211.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-50.5452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (188.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-58.8782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-108.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-68.5123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-250.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-37.9699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-39.5515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-152.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=COOC(=O)F(8762)'],
    products = ['CO2(14)', 'O=C(O)F(635)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(62.6195,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(13)', 'O=C(F)OO(2257)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.27e-07,'m^3/(mol*s)'), n=3.7, Ea=(159.122,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs_2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->H_1COCbCdCtHNOSSidSis->O',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs_2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->H_1COCbCdCtHNOSSidSis->O"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CO(13)', 'O=COOF(8964)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(80.8577,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C]1OO[CH]OO1(8965)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)OO[C]=O(8800)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C]OO[C](O)F(8760)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'O=[C]OOC=O(8966)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)F(459)', '[O]C=O(322)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.24e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC(=O)F(277)', 'HCO(15)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.5e+06,'m^3/(mol*s)'), n=-4.57487e-10, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_N-Sp-4R!H-2C_N-4R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_N-Sp-4R!H-2C_N-4R!H->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CFO(51)', '[O]OC=O(320)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'O=[C]OOC(=O)F(8799)'],
    products = ['O=COOC(=O)F(8762)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.31504e+12,'m^3/(mol*s)'), n=-2.20453, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O"""),
)

network(
    label = 'PDepNetwork #2485',
    isomers = [
        'O=COOC(=O)F(8762)',
    ],
    reactants = [
        ('CO2(14)', 'O=C(O)F(635)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2485',
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

