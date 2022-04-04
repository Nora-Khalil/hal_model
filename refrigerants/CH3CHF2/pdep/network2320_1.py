species(
    label = '[O]OC(F)=COOC(F)F(5911)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-579.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.325958,'amu*angstrom^2'), symmetry=1, barrier=(7.49441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325525,'amu*angstrom^2'), symmetry=1, barrier=(7.48447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325526,'amu*angstrom^2'), symmetry=1, barrier=(7.48447,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0773312,0.100144,-0.000183794,1.76058e-07,-6.44516e-11,-69595,34.8264], Tmin=(100,'K'), Tmax=(830.983,'K')), NASAPolynomial(coeffs=[6.91769,0.0386331,-2.11645e-05,4.19169e-09,-2.92519e-13,-69745,9.0303], Tmin=(830.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=COO(3391)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {2,S} {9,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {2,S} {7,D} {8,S}
7 C u0 p0 c0 {1,S} {3,S} {6,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-132.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.363982,'amu*angstrom^2'), symmetry=1, barrier=(8.36865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363921,'amu*angstrom^2'), symmetry=1, barrier=(8.36727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363788,'amu*angstrom^2'), symmetry=1, barrier=(8.3642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4277,0.0666494,-0.000124382,1.20591e-07,-4.42134e-11,-15895.7,25.1687], Tmin=(100,'K'), Tmax=(844.524,'K')), NASAPolynomial(coeffs=[5.19556,0.0267509,-1.43476e-05,2.81e-09,-1.94452e-13,-15745.7,12.2837], Tmin=(844.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
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
    label = '[O]OC(F)=COOF(3842)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {2,S} {3,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {3,S} {8,D} {9,S}
8 C u0 p0 c0 {1,S} {4,S} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-37.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,277,555,632,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180,1402.2],'cm^-1')),
        HinderedRotor(inertia=(0.275998,'amu*angstrom^2'), symmetry=1, barrier=(6.34574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275938,'amu*angstrom^2'), symmetry=1, barrier=(6.34436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275861,'amu*angstrom^2'), symmetry=1, barrier=(6.34258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06398,0.0763698,-0.000147358,1.43181e-07,-5.23298e-11,-4422.23,27.8552], Tmin=(100,'K'), Tmax=(845.402,'K')), NASAPolynomial(coeffs=[5.86452,0.0278442,-1.54602e-05,3.0532e-09,-2.11746e-13,-4311.51,10.957], Tmin=(845.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(C=O)OC(F)F(5985)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-962.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.279399,0.105698,-0.00018762,1.69777e-07,-5.91513e-11,-115654,31.6681], Tmin=(100,'K'), Tmax=(834.434,'K')), NASAPolynomial(coeffs=[10.709,0.0319608,-1.7207e-05,3.37669e-09,-2.34186e-13,-116754,-14.962], Tmin=(834.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-962.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFHO) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = 'O=C(F)[CH]OOC(F)F(5986)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
8  C u1 p0 c0 {5,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-814.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,519,593,830,1115,1166,1389,1413,3114,3025,407.5,1350,352.5,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.232195,'amu*angstrom^2'), symmetry=1, barrier=(5.33862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10253,'amu*angstrom^2'), symmetry=1, barrier=(48.3412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09901,'amu*angstrom^2'), symmetry=1, barrier=(48.2604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09916,'amu*angstrom^2'), symmetry=1, barrier=(48.2639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.593921,0.0820497,-0.000123464,1.01935e-07,-3.4312e-11,-97787.8,28.5174], Tmin=(100,'K'), Tmax=(723.536,'K')), NASAPolynomial(coeffs=[9.93481,0.0304128,-1.6419e-05,3.31024e-09,-2.36763e-13,-99139.6,-13.5248], Tmin=(723.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFHO) + group(COCsFO) + radical(OCJC=O)"""),
)

species(
    label = '[O]C(F)F(615)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-426.241,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(66.9995,'amu')),
        NonlinearRotor(inertia=([47.2408,48.1768,88.2639],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([480.226,502.831,615.399,942.972,1118.44,1159.86,1274.53,1324.03,2842.12],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07235,-0.00674902,8.67846e-05,-1.53143e-07,8.6196e-11,-51263.6,9.19392], Tmin=(10,'K'), Tmax=(577.889,'K')), NASAPolynomial(coeffs=[2.80054,0.016572,-1.1432e-05,3.6345e-09,-4.33761e-13,-51359,12.5348], Tmin=(577.889,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1=COOO1(3394)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u0 p0 c0 {2,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-63.0576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45648,-0.00296079,7.42328e-05,-9.36209e-08,3.49566e-11,-7551.03,17.4583], Tmin=(100,'K'), Tmax=(965.686,'K')), NASAPolynomial(coeffs=[9.89799,0.00855895,-2.99883e-06,6.60762e-10,-5.63624e-14,-10576.4,-22.616], Tmin=(965.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.0576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclopentane)"""),
)

species(
    label = 'F[C]1OOC1OOC(F)F(5987)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {5,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-619.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478233,0.0853492,-0.000126039,1.04418e-07,-3.56694e-11,-74393.1,32.2757], Tmin=(100,'K'), Tmax=(710.495,'K')), NASAPolynomial(coeffs=[9.498,0.0345696,-1.88342e-05,3.82728e-09,-2.75405e-13,-75674.8,-8.15666], Tmin=(710.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-619.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'FC(F)OO[CH]C1(F)OO1(5988)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {10,S}
8  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
10 C u1 p0 c0 {7,S} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-660.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.17462,0.092023,-0.000139026,1.14625e-07,-3.84231e-11,-79340.8,30.9338], Tmin=(100,'K'), Tmax=(727.372,'K')), NASAPolynomial(coeffs=[10.8479,0.0333144,-1.79287e-05,3.60819e-09,-2.57735e-13,-80893.1,-17.1588], Tmin=(727.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-660.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsFFHO) + ring(O2s-O2s-Cs(F)(C)) + radical(CCsJOOC)"""),
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
    label = '[O]OC(F)=COO[CH]F(3849)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {3,S} {8,D} {10,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u1 p0 c0 {2,S} {4,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-148.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,580,1155,1237,1373,3147,180],'cm^-1')),
        HinderedRotor(inertia=(0.311625,'amu*angstrom^2'), symmetry=1, barrier=(7.16488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311237,'amu*angstrom^2'), symmetry=1, barrier=(7.15595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311776,'amu*angstrom^2'), symmetry=1, barrier=(7.16835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397958,0.0939674,-0.000180147,1.76178e-07,-6.48909e-11,-17709.9,32.9606], Tmin=(100,'K'), Tmax=(842.788,'K')), NASAPolynomial(coeffs=[5.6855,0.0364657,-2.01287e-05,3.97579e-09,-2.76112e-13,-17450.2,15.1835], Tmin=(842.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFHHO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O]O[C]=COOC(F)F(5989)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-181.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.516513,'amu*angstrom^2'), symmetry=1, barrier=(11.8757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516601,'amu*angstrom^2'), symmetry=1, barrier=(11.8777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516535,'amu*angstrom^2'), symmetry=1, barrier=(11.8762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516587,'amu*angstrom^2'), symmetry=1, barrier=(11.8773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396406,0.0868362,-0.000145082,1.25238e-07,-4.2117e-11,-21713.2,34.798], Tmin=(100,'K'), Tmax=(830.137,'K')), NASAPolynomial(coeffs=[10.8678,0.0249974,-1.27761e-05,2.46868e-09,-1.70261e-13,-23059.6,-11.409], Tmin=(830.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]O[C](F)C=O(1584)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {2,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-189.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,280,501,1494,1531,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(4.16547,'amu*angstrom^2'), symmetry=1, barrier=(95.7722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.22148,'amu*angstrom^2'), symmetry=1, barrier=(97.0601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35167,0.0412352,-6.6895e-05,6.24069e-08,-2.28565e-11,-22685.8,17.8037], Tmin=(100,'K'), Tmax=(824.824,'K')), NASAPolynomial(coeffs=[4.73144,0.0198291,-1.00257e-05,1.94134e-09,-1.3458e-13,-22742.8,8.81534], Tmin=(824.824,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(ROOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'CHF2(81)',
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
    label = '[O]OC=C(F)O[O](3398)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u1 p2 c0 {2,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {2,S} {7,D} {8,S}
7 C u0 p0 c0 {1,S} {3,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (19.1511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([470,515,1100,1170,900,1100,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.233153,'amu*angstrom^2'), symmetry=1, barrier=(5.36065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234516,'amu*angstrom^2'), symmetry=1, barrier=(5.39199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68886,0.0630752,-0.000129038,1.29573e-07,-4.76594e-11,2374.86,24.991], Tmin=(100,'K'), Tmax=(871.403,'K')), NASAPolynomial(coeffs=[3.58325,0.025133,-1.33819e-05,2.57542e-09,-1.74843e-13,3155.11,22.4838], Tmin=(871.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.1511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(O2s-OsH) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)F(2933)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-427.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114],'cm^-1')),
        HinderedRotor(inertia=(0.586675,'amu*angstrom^2'), symmetry=1, barrier=(13.4888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3217.18,'J/mol'), sigma=(5.19785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.52 K, Pc=51.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92725,0.00470123,7.24333e-05,-1.66658e-07,1.13367e-10,-51470.2,11.559], Tmin=(10,'K'), Tmax=(497.316,'K')), NASAPolynomial(coeffs=[3.77619,0.0196241,-1.39226e-05,4.52949e-09,-5.50671e-13,-51624.7,10.478], Tmin=(497.316,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-427.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=C(F)O[O](1081)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u1 p0 c0 {4,D} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (157.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,293,496,537,1218,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.441569,'amu*angstrom^2'), symmetry=1, barrier=(10.1525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3364.7,'J/mol'), sigma=(5.40718,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.56 K, Pc=48.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00596,0.051545,-0.000104345,9.71572e-08,-3.31603e-11,19053.1,15.8311], Tmin=(100,'K'), Tmax=(897.319,'K')), NASAPolynomial(coeffs=[7.11388,0.0105561,-5.37028e-06,9.95721e-10,-6.48097e-14,18869.9,-4.17116], Tmin=(897.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Cdj(Cd-F1sO2s)(H))"""),
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
    label = 'F[C]=COOC(F)F(5990)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
7  C u0 p0 c0 {5,S} {8,D} {10,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-395.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,167,640,1190,180],'cm^-1')),
        HinderedRotor(inertia=(0.517656,'amu*angstrom^2'), symmetry=1, barrier=(11.9019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20669,'amu*angstrom^2'), symmetry=1, barrier=(50.736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.517073,'amu*angstrom^2'), symmetry=1, barrier=(11.8885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.871402,0.0768913,-0.000130579,1.18415e-07,-4.21325e-11,-47407.3,27.1259], Tmin=(100,'K'), Tmax=(806.483,'K')), NASAPolynomial(coeffs=[8.23687,0.0275003,-1.47972e-05,2.93438e-09,-2.0597e-13,-48177.1,-4.23135], Tmin=(806.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-395.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(CsFFHO) + group(Cds-CdsOsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Cdj(Cd-O2sH)(F1s))"""),
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
    label = '[O]OC(F)=COO[C](F)F(5991)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {1,S} {6,S} {8,D}
10 C u1 p0 c0 {2,S} {3,S} {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-361.695,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.280203,'amu*angstrom^2'), symmetry=1, barrier=(6.44241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28059,'amu*angstrom^2'), symmetry=1, barrier=(6.45132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280804,'amu*angstrom^2'), symmetry=1, barrier=(6.45624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.59688,'amu*angstrom^2'), symmetry=1, barrier=(59.7074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.139474,0.0977287,-0.000180032,1.69713e-07,-6.11113e-11,-43375.3,35.1018], Tmin=(100,'K'), Tmax=(834.245,'K')), NASAPolynomial(coeffs=[8.10731,0.0337434,-1.86285e-05,3.68716e-09,-2.5681e-13,-43807.5,3.48262], Tmin=(834.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-361.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = '[O]OC(F)=[C]OOC(F)F(5992)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-339.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,1685,370,180,1709.2],'cm^-1')),
        HinderedRotor(inertia=(0.282711,'amu*angstrom^2'), symmetry=1, barrier=(6.50007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282933,'amu*angstrom^2'), symmetry=1, barrier=(6.50519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282511,'amu*angstrom^2'), symmetry=1, barrier=(6.49548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281884,'amu*angstrom^2'), symmetry=1, barrier=(6.48108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270087,0.102107,-0.000211395,2.15721e-07,-8.0992e-11,-40773.3,37.2627], Tmin=(100,'K'), Tmax=(854.549,'K')), NASAPolynomial(coeffs=[2.8954,0.0425521,-2.38894e-05,4.71413e-09,-3.25785e-13,-39496.2,35.1076], Tmin=(854.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-339.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CJO)"""),
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
    label = '[O]OC#COOC(F)F(5993)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
8  C u0 p0 c0 {4,S} {9,T}
9  C u0 p0 c0 {5,S} {8,T}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-130.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114,2100,2250,500,550,180],'cm^-1')),
        HinderedRotor(inertia=(1.45551,'amu*angstrom^2'), symmetry=1, barrier=(33.465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46971,'amu*angstrom^2'), symmetry=1, barrier=(33.7915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05141,'amu*angstrom^2'), symmetry=1, barrier=(24.1741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43879,'amu*angstrom^2'), symmetry=1, barrier=(33.0807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387015,0.0879023,-0.000151357,1.34243e-07,-4.65445e-11,-15557.8,30.9166], Tmin=(100,'K'), Tmax=(806.533,'K')), NASAPolynomial(coeffs=[10.4902,0.026008,-1.43219e-05,2.85074e-09,-2.00109e-13,-16804.1,-13.2762], Tmin=(806.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(O2s-OsCt) + group(O2s-OsH) + group(CsFFHO) + group(Ct-CtOs) + group(Ct-CtOs) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)=[C]OOC(F)F(5994)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {5,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-491.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,519,593,830,1115,1166,1389,1413,3114,293,496,537,1218,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.168762,'amu*angstrom^2'), symmetry=1, barrier=(3.88017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168738,'amu*angstrom^2'), symmetry=1, barrier=(3.87962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168961,'amu*angstrom^2'), symmetry=1, barrier=(3.88475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168984,'amu*angstrom^2'), symmetry=1, barrier=(3.88527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.29336,'amu*angstrom^2'), symmetry=1, barrier=(52.7288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00949041,0.105683,-0.000206793,2.06913e-07,-7.76861e-11,-59043.9,37.438], Tmin=(100,'K'), Tmax=(837.661,'K')), NASAPolynomial(coeffs=[4.47791,0.0442216,-2.48853e-05,4.95593e-09,-3.45996e-13,-58384.8,25.0745], Tmin=(837.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-491.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'OOC(F)=COO[C](F)F(5995)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {1,S} {6,S} {8,D}
10 C u1 p0 c0 {2,S} {3,S} {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-513.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.776469,'amu*angstrom^2'), symmetry=1, barrier=(17.8525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.38419,'amu*angstrom^2'), symmetry=1, barrier=(54.8173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205267,'amu*angstrom^2'), symmetry=1, barrier=(4.71949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20528,'amu*angstrom^2'), symmetry=1, barrier=(4.71979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.38824,'amu*angstrom^2'), symmetry=1, barrier=(54.9104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0986001,0.101031,-0.000174461,1.59698e-07,-5.73554e-11,-61646.8,35.197], Tmin=(100,'K'), Tmax=(798.606,'K')), NASAPolynomial(coeffs=[9.5771,0.035613,-1.97434e-05,3.9577e-09,-2.79446e-13,-62651.5,-5.92168], Tmin=(798.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-513.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFHO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'F[CH]OOC=C(F)OOF(5996)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {3,S} {6,S}
8  C u0 p0 c0 {4,S} {9,D} {11,S}
9  C u0 p0 c0 {1,S} {6,S} {8,D}
10 C u1 p0 c0 {2,S} {5,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-204.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,3010,987.5,1337.5,450,1655,326,540,652,719,1357,580,1155,1237,1373,3147],'cm^-1')),
        HinderedRotor(inertia=(1.40441,'amu*angstrom^2'), symmetry=1, barrier=(32.2902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990812,'amu*angstrom^2'), symmetry=1, barrier=(22.7807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.897095,'amu*angstrom^2'), symmetry=1, barrier=(20.626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.218294,0.107168,-0.000198205,1.89608e-07,-6.95987e-11,-24507.3,35.7935], Tmin=(100,'K'), Tmax=(819.972,'K')), NASAPolynomial(coeffs=[7.88705,0.0393166,-2.22893e-05,4.47334e-09,-3.14674e-13,-24884.8,4.10249], Tmin=(819.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFHHO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(CsF1sHO2s)"""),
)

species(
    label = 'FOO[C]=COOC(F)F(5997)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {3,S} {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-238.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,277,555,632,519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.04838,'amu*angstrom^2'), symmetry=1, barrier=(24.1044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0489,'amu*angstrom^2'), symmetry=1, barrier=(24.1163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04728,'amu*angstrom^2'), symmetry=1, barrier=(24.079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0778533,0.0957955,-0.000143979,1.05196e-07,-2.70481e-11,-28523.4,36.6044], Tmin=(100,'K'), Tmax=(625.473,'K')), NASAPolynomial(coeffs=[12.7511,0.0284445,-1.53046e-05,3.05746e-09,-2.16667e-13,-30376.7,-20.7314], Tmin=(625.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFFHO) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
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
    E0 = (-319.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-28.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (351.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-256.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-332.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-250.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-214.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-302.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (163.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (130.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-329.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1.10965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-31.4053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-164.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (88.7784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (110.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (113.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-208.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-212.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (94.6821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (93.2289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)=COOC(F)F(5911)'],
    products = ['CF2O(49)', '[O]OC(F)C=O(1537)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(21.4708,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(43)', '[O]OC(F)=COO(3391)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(69.7321,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', '[O]OC(F)=COOF(3842)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(12.0472,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC(F)=COOC(F)F(5911)'],
    products = ['[O]OC(F)(C=O)OC(F)F(5985)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(85.022,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C(F)[CH]OOC(F)F(5986)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)=COOC(F)F(5911)'],
    products = ['[O]C(F)F(615)', 'FC1=COOO1(3394)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(90.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 86.6 to 90.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OC(F)=COOC(F)F(5911)'],
    products = ['F[C]1OOC1OOC(F)F(5987)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC(F)=COOC(F)F(5911)'],
    products = ['FC(F)OO[CH]C1(F)OO1(5988)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(38.4928,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[O]OC(F)=COO[CH]F(3849)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]O[C]=COOC(F)F(5989)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)F(615)', '[O]O[C](F)C=O(1584)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(47.1603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHF2(81)', '[O]OC=C(F)O[O](3398)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC(F)F(2933)', '[CH]=C(F)O[O](1081)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.14976e+13,'m^3/(mol*s)'), n=-2.27223, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O2(2)', 'F[C]=COOC(F)F(5990)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[O]OC(F)=COO[C](F)F(5991)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]OC(F)=[C]OOC(F)F(5992)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', '[O]OC#COOC(F)F(5993)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(286.17,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction18',
    reactants = ['OOC(F)=[C]OOC(F)F(5994)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OOC(F)=COO[C](F)F(5995)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2000,'s^-1'), n=1.9, Ea=(62.3416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H_OOCs4;C_rad_out_single;XH_out] for rate rule [R7H_OOCs4;C_rad_out_noH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]OOC=C(F)OOF(5996)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(60.9255,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction21',
    reactants = ['FOO[C]=COOC(F)F(5997)'],
    products = ['[O]OC(F)=COOC(F)F(5911)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(92.8166,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #2320',
    isomers = [
        '[O]OC(F)=COOC(F)F(5911)',
    ],
    reactants = [
        ('CF2O(49)', '[O]OC(F)C=O(1537)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #2320',
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

