species(
    label = 'O=C(OF)C1OC1F(2819)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-532.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2950,1000,261,493,600,1152,1365,1422,3097,180,1043.27,1043.27,1043.27,1043.27,1043.27,1043.27,1043.27,1043.27,1043.27,2266.38],'cm^-1')),
        HinderedRotor(inertia=(0.0378523,'amu*angstrom^2'), symmetry=1, barrier=(0.870299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0378523,'amu*angstrom^2'), symmetry=1, barrier=(0.870299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69142,0.0384278,-1.10698e-06,-2.86836e-08,1.36667e-11,-63948,24.4138], Tmin=(100,'K'), Tmax=(1026.8,'K')), NASAPolynomial(coeffs=[15.0788,0.0142227,-6.57293e-06,1.37213e-09,-1.0483e-13,-68170.5,-47.7006], Tmin=(1026.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-532.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cs-Cs(F)-O2s)"""),
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
    label = 'O=CC(=O)OF(1291)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,D} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-347.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2782.5,750,1395,475,1775,1000,443.665,444.054,444.153,444.627,444.964],'cm^-1')),
        HinderedRotor(inertia=(0.332923,'amu*angstrom^2'), symmetry=1, barrier=(46.6435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331199,'amu*angstrom^2'), symmetry=1, barrier=(46.6382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.65562,0.0341149,-3.02436e-05,1.3375e-08,-2.52359e-12,-41726.2,14.5066], Tmin=(100,'K'), Tmax=(1189.66,'K')), NASAPolynomial(coeffs=[7.40424,0.0181487,-1.01124e-05,2.09375e-09,-1.52915e-13,-42856,-9.22741], Tmin=(1189.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = 'O-2(136)',
    structure = adjacencyList("""1 O u0 p3 c0
"""),
    E0 = (432.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-1.84483e-15,2.71425e-18,-1.30028e-21,1.91033e-25,51997.4,2.99252], Tmin=(100,'K'), Tmax=(3598.68,'K')), NASAPolynomial(coeffs=[2.5,-2.82485e-12,1.07037e-15,-1.78888e-19,1.11248e-23,51997.4,2.99252], Tmin=(3598.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.331,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(C=CF)OF(666)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {6,S} {7,D} {8,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u0 p0 c0 {1,S} {5,D} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-344.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3010,987.5,1337.5,450,1655,194,682,905,1196,1383,3221,180,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.294089,'amu*angstrom^2'), symmetry=1, barrier=(6.76169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292402,'amu*angstrom^2'), symmetry=1, barrier=(6.7229,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21694,0.0441445,-7.46461e-05,7.01091e-08,-2.62074e-11,-41419.2,17.5915], Tmin=(100,'K'), Tmax=(766.145,'K')), NASAPolynomial(coeffs=[5.89325,0.018029,-9.96401e-06,2.03336e-09,-1.45831e-13,-41779.3,2.16052], Tmin=(766.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFH)"""),
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
    label = 'FOC1OC1F(3035)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {2,S} {5,S}
5 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-337.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,261,493,600,1152,1365,1422,3097,605.066,605.076,605.077,605.077,605.077,605.081],'cm^-1')),
        HinderedRotor(inertia=(0.000460447,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2875,0.031236,-1.44578e-05,-3.30169e-09,3.01636e-12,-40501.8,15.7067], Tmin=(100,'K'), Tmax=(1125.31,'K')), NASAPolynomial(coeffs=[10.7677,0.0130543,-6.16728e-06,1.23308e-09,-8.9691e-14,-43167.7,-29.572], Tmin=(1125.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(Cs(O2)-O2s-Cs(F))"""),
)

species(
    label = 'OF(147)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=C1OC1F(3036)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5 C u0 p0 c0 {2,S} {4,S} {6,D}
6 C u0 p0 c0 {3,D} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-235.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([236,527,855,1015,1182,1348,3236,2120,512.5,787.5,206.545,207.468,208.118,208.404,209.133],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0371,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88502,0.0333929,-2.2333e-06,-3.73507e-08,2.12717e-11,-28182.3,11.4508], Tmin=(100,'K'), Tmax=(916.174,'K')), NASAPolynomial(coeffs=[18.6553,-0.00520401,4.27501e-06,-8.39504e-10,5.35438e-14,-32708.3,-75.9177], Tmin=(916.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-235.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFHO) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + ring(methyleneoxirane)"""),
)

species(
    label = 'FOC1=CC(F)OO1(3037)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {2,S} {8,S}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
7  C u0 p0 c0 {6,S} {8,D} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-226.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,236,527,855,1015,1182,1348,3236,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84115,0.0534165,-5.67698e-05,3.35431e-08,-8.54361e-12,-27118.4,19.9532], Tmin=(100,'K'), Tmax=(916.271,'K')), NASAPolynomial(coeffs=[7.70626,0.0278118,-1.48523e-05,3.04391e-09,-2.21902e-13,-28193.2,-7.82958], Tmin=(916.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(12dioxolene)"""),
)

species(
    label = 'FOC1=COC(F)O1(3038)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {5,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {2,S} {7,S}
6  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-518.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,509,613,660,1171,1360,1414,3084,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.994448,0.0413511,3.12726e-05,-9.47016e-08,4.59356e-11,-62282.8,18.2176], Tmin=(100,'K'), Tmax=(921.86,'K')), NASAPolynomial(coeffs=[27.174,-0.00792027,6.78124e-06,-1.30024e-09,7.96672e-14,-69842.7,-120.778], Tmin=(921.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-518.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsFHOO) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclopentane)"""),
)

species(
    label = 'O=C([CH]O[CH]F)OF(3039)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {7,D}
6  C u1 p0 c0 {3,S} {7,S} {9,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u1 p0 c0 {1,S} {3,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-268.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3025,407.5,1350,352.5,580,1155,1237,1373,3147,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825315,0.0753751,-0.00010858,8.02562e-08,-2.36471e-11,-32218.4,25.5594], Tmin=(100,'K'), Tmax=(829.606,'K')), NASAPolynomial(coeffs=[11.9462,0.0217542,-1.16273e-05,2.34421e-09,-1.68158e-13,-34063.5,-26.0147], Tmin=(829.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(CsFHHO) + group(Cds-OdCsOs) + radical(CCsJOCs) + radical(Csj(F1s)(O2s-Cs)(H))"""),
)

species(
    label = '[O]C([CH]F)C(=O)OF(3040)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {5,D} {6,S}
8  C u1 p0 c0 {1,S} {6,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-254.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,334,575,1197,1424,3202,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36148,0.062508,-7.51297e-05,4.71026e-08,-1.20197e-11,-30486.7,27.5675], Tmin=(100,'K'), Tmax=(943.227,'K')), NASAPolynomial(coeffs=[10.7375,0.022747,-1.18991e-05,2.41214e-09,-1.74731e-13,-32255.4,-17.1183], Tmin=(943.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFHH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(CsCsF1sH)"""),
)

species(
    label = '[O]C(F)[CH]C(=O)OF(3041)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {3,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u1 p0 c0 {6,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {5,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-293.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,391,562,707,872,1109,1210,1289,3137,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45957,0.0612584,-7.95428e-05,5.80679e-08,-1.76348e-11,-35220.8,26.2779], Tmin=(100,'K'), Tmax=(792.993,'K')), NASAPolynomial(coeffs=[8.35322,0.0264839,-1.37612e-05,2.7628e-09,-1.98349e-13,-36314,-5.38076], Tmin=(792.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-293.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsHH) + group(CsCFHO) + group(Cds-OdCsOs) + radical(O2sj(Cs-CsF1sH)) + radical(CCJCO)"""),
)

species(
    label = '[O]C(OF)[C]1OC1F(3042)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-135.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,487,638,688,1119,1325,1387,3149,1380,1390,370,380,2900,435,180,180,180,1313.87,1314.03,1314.06,1314.71],'cm^-1')),
        HinderedRotor(inertia=(0.189427,'amu*angstrom^2'), symmetry=1, barrier=(4.35529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189134,'amu*angstrom^2'), symmetry=1, barrier=(4.34856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4349,0.065379,-0.000111218,1.08439e-07,-4.12717e-11,-16181.9,27.4435], Tmin=(100,'K'), Tmax=(809.125,'K')), NASAPolynomial(coeffs=[4.27271,0.0331249,-1.76373e-05,3.49679e-09,-2.4579e-13,-16044.6,18.0403], Tmin=(809.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(O2s-Cs(F)-Cs(C)) + radical(CCOJ) + radical(C2CsJO)"""),
)

species(
    label = '[O]C(OF)C1O[C]1F(3043)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
8  C u1 p0 c0 {1,S} {3,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-92.2491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2950,1000,1380,1390,370,380,2900,435,395,473,707,1436,305.682,305.682,305.683,864.464,1975.6,1975.6,1975.6,1975.6],'cm^-1')),
        HinderedRotor(inertia=(0.521568,'amu*angstrom^2'), symmetry=1, barrier=(34.5843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.521568,'amu*angstrom^2'), symmetry=1, barrier=(34.5843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70392,0.0572019,-6.94513e-05,4.31976e-08,-7.71038e-12,-11018.7,26.1937], Tmin=(100,'K'), Tmax=(571.218,'K')), NASAPolynomial(coeffs=[6.78301,0.0291317,-1.54251e-05,3.11869e-09,-2.24562e-13,-11721.2,3.46374], Tmin=(571.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.2491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(O2s-Cs(F)-Cs(C)) + radical(CCOJ) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'O[C](OF)C1O[C]1F(3044)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {8,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {7,S} {10,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u1 p0 c0 {1,S} {3,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-112.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,2950,1000,360,370,350,395,473,707,1436,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34544,0.0615966,-7.40249e-05,4.55059e-08,-1.12333e-11,-13462.7,27.2793], Tmin=(100,'K'), Tmax=(980.406,'K')), NASAPolynomial(coeffs=[11.5879,0.0198073,-1.00872e-05,2.02811e-09,-1.46436e-13,-15471,-21.9318], Tmin=(980.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(CsCFHO) + ring(O2s-Cs(F)-Cs(C)) + radical(Cs_P) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'OC(OF)=C1OC1F(3045)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {2,S} {8,S}
5  O u0 p2 c0 {8,S} {10,S}
6  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {5,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-367.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,236,527,855,1015,1182,1348,3236,350,440,435,1725,376.886,376.98,376.985,377.014,377.291,1325.88],'cm^-1')),
        HinderedRotor(inertia=(0.00118716,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117177,'amu*angstrom^2'), symmetry=1, barrier=(11.8141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968617,0.0692764,-9.17566e-05,6.06048e-08,-1.57285e-11,-44041.7,22.9928], Tmin=(100,'K'), Tmax=(944.564,'K')), NASAPolynomial(coeffs=[13.3982,0.016641,-8.17115e-06,1.61189e-09,-1.15017e-13,-46389.9,-36.2642], Tmin=(944.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + ring(Cs(F)-Cd(Cd)-O2s)"""),
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
    label = 'O=C(OF)C1[CH]O1(3046)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {5,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-131.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,3150,900,1100,180,1065.51,1065.51,1065.51,1065.51,1065.51,1065.51,1065.51,1065.51,1065.51,1065.51,1065.51,2254.91],'cm^-1')),
        HinderedRotor(inertia=(0.0526439,'amu*angstrom^2'), symmetry=1, barrier=(1.21039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0526439,'amu*angstrom^2'), symmetry=1, barrier=(1.21039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72165,0.0374291,4.92312e-06,-4.79089e-08,2.61438e-11,-15673.3,21.1086], Tmin=(100,'K'), Tmax=(886.144,'K')), NASAPolynomial(coeffs=[17.0939,0.00354777,2.16984e-06,-6.19086e-10,4.52737e-14,-19791.8,-59.0615], Tmin=(886.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Ethylene_oxide) + radical(CCsJO)"""),
)

species(
    label = '[O]C(=O)C1OC1F(3047)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-441.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,261,493,600,1152,1365,1422,3097,600,948.753,948.753,948.753,948.753,948.753,948.753,948.753,948.753,948.753,2337.86],'cm^-1')),
        HinderedRotor(inertia=(0.00243889,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55098,0.0229986,1.4796e-05,-3.10175e-08,1.16244e-11,-53013.8,21.3295], Tmin=(100,'K'), Tmax=(1085.71,'K')), NASAPolynomial(coeffs=[9.52391,0.018781,-9.04265e-06,1.83631e-09,-1.3513e-13,-55793.4,-18.7124], Tmin=(1085.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cs-Cs(F)-O2s) + radical(CCOJ)"""),
)

species(
    label = '[O]F(125)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2889.11,'J/mol'), sigma=(4.75593,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=451.27 K, Pc=60.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]C1OC1F(3048)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-241.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,261,493,600,1152,1365,1422,3097,1855,455,950,621.899,622.677,623.715,623.977,2742.72],'cm^-1')),
        HinderedRotor(inertia=(0.00043458,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0451,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63427,0.0183071,2.51201e-05,-5.01088e-08,2.1412e-11,-28966.2,20.3661], Tmin=(100,'K'), Tmax=(959.132,'K')), NASAPolynomial(coeffs=[12.5358,0.00625825,-1.77332e-06,3.74493e-10,-3.28153e-14,-32210.7,-34.0025], Tmin=(959.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsH) + ring(Cs-Cs(F)-O2s) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]OF(146)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {2,S} {3,D}
"""),
    E0 = (-35.6837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.77736,'amu*angstrom^2'), symmetry=1, barrier=(40.8651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98172,0.0254564,-4.70449e-05,4.25373e-08,-1.44425e-11,-4258,9.53187], Tmin=(100,'K'), Tmax=(878.476,'K')), NASAPolynomial(coeffs=[5.67747,0.00648674,-3.22257e-06,6.05546e-10,-4.04972e-14,-4473.3,-1.65391], Tmin=(878.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.6837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical((O)CJOC)"""),
)

species(
    label = 'FC1[CH]O1(2761)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-64.5209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,2950,1000,592.694,1288.9,1288.9],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09538,-0.00819916,9.64624e-05,-1.59005e-07,8.41279e-11,-7758.89,9.44595], Tmin=(10,'K'), Tmax=(607.315,'K')), NASAPolynomial(coeffs=[2.06726,0.0207746,-1.36692e-05,4.22801e-09,-4.95559e-13,-7800.53,15.8481], Tmin=(607.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-64.5209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC1[CH]O1""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(OF)[C]1OC1F(3049)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {2,S} {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
7 C u1 p0 c0 {3,S} {6,S} {8,S}
8 C u0 p0 c0 {4,S} {5,D} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-352.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,487,638,688,1119,1325,1387,3149,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6191,0.0450609,-3.63824e-05,1.35686e-08,-1.97313e-12,-42288.6,25.8055], Tmin=(100,'K'), Tmax=(1643.36,'K')), NASAPolynomial(coeffs=[15.4435,0.0114118,-5.66866e-06,1.10886e-09,-7.76664e-14,-46832.3,-47.7565], Tmin=(1643.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-352.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cs-Cs(F)-O2s) + radical(C2CsJO)"""),
)

species(
    label = 'O=C(OF)C1O[C]1F(3050)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {2,S} {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {1,S} {3,S} {6,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-337.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2950,1000,395,473,707,1436,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56089,0.0467122,-3.96796e-05,1.57282e-08,-2.43335e-12,-40549,26.3906], Tmin=(100,'K'), Tmax=(1551.74,'K')), NASAPolynomial(coeffs=[15.2143,0.0115173,-5.65829e-06,1.11187e-09,-7.8525e-14,-44786.3,-45.4783], Tmin=(1551.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=C(OF)C1=CO1(3034)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {6,D} {7,S}
6 C u0 p0 c0 {2,S} {5,D} {8,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-3.53532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57768,0.0557966,-7.52496e-05,4.96174e-08,-1.28241e-11,-340.056,19.2663], Tmin=(100,'K'), Tmax=(948.361,'K')), NASAPolynomial(coeffs=[11.8726,0.0123748,-6.57054e-06,1.3383e-09,-9.71578e-14,-2292.72,-29.8549], Tmin=(948.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.53532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2sCF) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cyclopropene)"""),
)

species(
    label = 'F2(77)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C1=CO1(3051)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u0 p0 c0 {1,S} {2,D} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (259.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31417,0.0119135,-5.00127e-07,-9.18336e-09,4.80348e-12,31257.7,8.09654], Tmin=(100,'K'), Tmax=(945.068,'K')), NASAPolynomial(coeffs=[6.90182,0.00459354,-1.36476e-06,2.32146e-10,-1.65655e-14,30228.3,-10.8673], Tmin=(945.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""oxirene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    E0 = (-308.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-12.2599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (245.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-94.7284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-1.16848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (186.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (25.5639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-102.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-92.8645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-127.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (45.5074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (129.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (54.1063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-62.5884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (100.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-209.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (19.7964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (58.2415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (17.8806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (32.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (28.5825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (32.4762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(OF)C1OC1F(2819)'],
    products = ['HF(38)', 'CO2(14)', 'FC1=CO1(1231)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(65.7735,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', 'O=CC(=O)OF(1291)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.81575e-05,'m^3/(mol*s)'), n=2.76922, Ea=(37.8393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;mb_carbonyl] for rate rule [carbene;mb_carbonyl_HDe]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O-2(136)', 'O=C(C=CF)OF(666)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.54e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [o_atom_singlet;mb_db_twocdisub]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CO(13)', 'FOC1OC1F(3035)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(202.875,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(147)', 'O=C=C1OC1F(3036)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1450,'cm^3/(mol*s)'), n=2.8, Ea=(170.707,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 6 used for cco_Nd2;H_OH
Exact match found for rate rule [cco_Nd2;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FOC1=CC(F)OO1(3037)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(254.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOC1=COC(F)O1(3038)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.09884e+62,'s^-1'), n=-13.6281, Ea=(386.057,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.1335689029548193, var=267.23618516426137, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_1R!H-inRing',), comment="""Estimated from node Root_1R!H-inRing"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C([CH]O[CH]F)OF(3039)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] for rate rule [R3_SS;C_rad_out_H/OneDe;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C([CH]F)C(=O)OF(3040)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_1H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)[CH]C(=O)OF(3041)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(OF)[C]1OC1F(3042)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(OF)C1O[C]1F(3043)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O[C](OF)C1O[C]1F(3044)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OC(OF)=C1OC1F(3045)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(146.039,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'O=C(OF)C1[CH]O1(3046)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', '[O]C(=O)C1OC1F(3047)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.11128e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]F(125)', 'O=[C]C1OC1F(3048)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]OF(146)', 'FC1[CH]O1(2761)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.17885e+10,'m^3/(mol*s)'), n=-0.943326, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing_Ext-2R-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing_Ext-2R-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'O=C(OF)[C]1OC1F(3049)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(176351,'m^3/(mol*s)'), n=-0.549379, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Ext-2BrCClFHNO-R_Ext-5R!H-R_Ext-5R!H-R_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'O=C(OF)C1O[C]1F(3050)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.37186e+06,'m^3/(mol*s)'), n=0.39882, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'O=C(OF)C1=CO1(3034)'],
    products = ['O=C(OF)C1OC1F(2819)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(154.785,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C(OF)C1OC1F(2819)'],
    products = ['F2(77)', 'CO2(14)', 'C1=CO1(3051)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(406.505,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

network(
    label = 'PDepNetwork #1269',
    isomers = [
        'O=C(OF)C1OC1F(2819)',
    ],
    reactants = [
        ('HF(38)', 'CO2(14)', 'FC1=CO1(1231)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1269',
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

