species(
    label = 'F[C](F)OC1[C](F)C=C1F(13400)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
7  C u1 p0 c0 {1,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {12,S}
10 C u1 p0 c0 {3,S} {4,S} {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-521.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,493,600,700,1144,1293,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856544,0.0697278,-6.83454e-05,3.30855e-08,-6.44207e-12,-62608.3,29.5558], Tmin=(100,'K'), Tmax=(1222.28,'K')), NASAPolynomial(coeffs=[14.8386,0.0239709,-1.21921e-05,2.45796e-09,-1.77685e-13,-66026.3,-40.7058], Tmin=(1222.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-521.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CsFFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
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
    label = 'FC1=CC(F)=C1(307)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {4,S} {5,D} {7,S}
4 C u0 p0 c0 {1,S} {3,S} {6,D}
5 C u0 p0 c0 {2,S} {3,D} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {8,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (29.0416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,217,343,449,587,660,812,790,914,798,948,180,1170.89,2106.02,2106.02],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3241.95,'J/mol'), sigma=(5.01646,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.38 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89011,0.00686793,0.000107406,-2.40335e-07,1.58082e-10,3500.24,9.34061], Tmin=(10,'K'), Tmax=(520.912,'K')), NASAPolynomial(coeffs=[4.12106,0.0279909,-1.93509e-05,6.26873e-09,-7.66374e-13,3165.53,5.39524], Tmin=(520.912,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(29.0416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCC(F)DC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)O[CH]C1(F)C=C1F(13438)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  C u1 p0 c0 {5,S} {6,S} {12,S}
10 C u1 p0 c0 {3,S} {4,S} {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-353.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,323,467,575,827,1418,2950,1000,3025,407.5,1350,352.5,493,600,700,1144,1293,265.611,265.616,265.617,1032.57,1032.58],'cm^-1')),
        HinderedRotor(inertia=(0.163579,'amu*angstrom^2'), symmetry=1, barrier=(8.18946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163567,'amu*angstrom^2'), symmetry=1, barrier=(8.18934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989143,'amu*angstrom^2'), symmetry=1, barrier=(49.5203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0796622,0.0936531,-0.000140984,1.07892e-07,-3.17549e-11,-42373.2,31.4719], Tmin=(100,'K'), Tmax=(684.747,'K')), NASAPolynomial(coeffs=[13.0404,0.0267348,-1.36555e-05,2.67767e-09,-1.88014e-13,-44354.3,-27.6531], Tmin=(684.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cs-CsOsHH) + group(CsFFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cd-Cs(F)(C)-Cd) + radical(CCsJOCs) + radical(Csj(F1s)(F1s)(O2s-Cs)) + longDistanceInteraction_cyclic(Cs(F)-Cds(F))"""),
)

species(
    label = 'F[C]F(156)',
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
    label = '[O]C1[C](F)C=C1F(1459)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
5 C u1 p0 c0 {1,S} {4,S} {7,S}
6 C u0 p0 c0 {2,S} {4,S} {7,D}
7 C u0 p0 c0 {5,S} {6,D} {9,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-76.9062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,323,467,575,827,1418,180,864.045,864.046,864.046,864.048,864.055,864.056,864.058],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92628,0.0338811,1.66589e-06,-3.10557e-08,1.50484e-11,-9164.58,20.4671], Tmin=(100,'K'), Tmax=(991.386,'K')), NASAPolynomial(coeffs=[14.6396,0.0101379,-4.09666e-06,8.52018e-10,-6.69499e-14,-13039.3,-47.5857], Tmin=(991.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.9062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CC(C)OJ) + radical(CsCdCsF1s)"""),
)

species(
    label = 'FC1=CC2(F)C1OC2(F)F(13401)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-740.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959239,0.0669916,-6.0029e-05,2.56858e-08,-4.41094e-12,-88988,20.5882], Tmin=(100,'K'), Tmax=(1367.7,'K')), NASAPolynomial(coeffs=[15.7422,0.0237576,-1.26134e-05,2.57392e-09,-1.86396e-13,-93031.8,-55.3599], Tmin=(1367.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-740.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCsCdF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + polycyclic(s2_4_4_ene_1)"""),
)

species(
    label = 'FC1=CC(F)=C1OC(F)F(13439)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {9,S}
8  C u0 p0 c0 {3,S} {7,D} {10,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-584.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441741,0.0847234,-0.000112045,7.70732e-08,-2.13833e-11,-70145.5,25.1457], Tmin=(100,'K'), Tmax=(874.704,'K')), NASAPolynomial(coeffs=[12.8311,0.028067,-1.48866e-05,3.02211e-09,-2.18562e-13,-72312.9,-32.9671], Tmin=(874.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-584.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFFHO) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cd-Cd-Cd-Cd(F))"""),
)

species(
    label = 'FC1=C=C(F)C1OC(F)F(13403)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u0 p0 c0 {4,S} {6,S} {10,D}
10 C u0 p0 c0 {8,D} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-478.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.496934,0.085906,-0.000139061,1.27274e-07,-4.69153e-11,-57447.5,28.6805], Tmin=(100,'K'), Tmax=(764.523,'K')), NASAPolynomial(coeffs=[7.79087,0.0364472,-1.98587e-05,4.00189e-09,-2.8506e-13,-58232.7,-2.3909], Tmin=(764.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-478.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(CsFFHO) + group(CdCddCF) + group(CdCddCF) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'F[C](F)OC1C2(F)[CH]C12F(13440)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
9  C u1 p0 c0 {6,S} {7,S} {12,S}
10 C u1 p0 c0 {3,S} {4,S} {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-392.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341932,0.0908944,-0.000155134,1.44778e-07,-5.31789e-11,-47067.2,30.2083], Tmin=(100,'K'), Tmax=(795.7,'K')), NASAPolynomial(coeffs=[7.64349,0.0367179,-2.00683e-05,4.02031e-09,-2.84262e-13,-47676,0.126641], Tmin=(795.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-392.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsFFHO) + polycyclic(s2_3_3_ane) + radical(bicyclo[1.1.0]butane-secondary) + radical(Csj(F1s)(F1s)(O2s-Cs)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C]1C2OC(F)(F)C1[C]2F(13441)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {8,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {9,S} {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
9  C u1 p0 c0 {3,S} {6,S} {7,S}
10 C u1 p0 c0 {4,S} {6,S} {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-549.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56536,0.057798,-4.3606e-05,1.52226e-08,-2.15788e-12,-65998.4,21.4937], Tmin=(100,'K'), Tmax=(1577.46,'K')), NASAPolynomial(coeffs=[13.6377,0.027186,-1.44974e-05,2.9208e-09,-2.08273e-13,-69807.1,-42.2515], Tmin=(1577.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-549.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFH) + group(CsCFFO) + polycyclic(s3_4_5_ane) + radical(CsCsCsF1s) + radical(CsCsCsF1s)"""),
)

species(
    label = 'F[C]1[CH]C2(F)C1OC2(F)F(13442)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
9  C u1 p0 c0 {6,S} {10,S} {12,S}
10 C u1 p0 c0 {4,S} {7,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-480.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.027,0.0731367,-9.67975e-05,7.94559e-08,-2.84566e-11,-57714,24.2868], Tmin=(100,'K'), Tmax=(666.181,'K')), NASAPolynomial(coeffs=[6.56413,0.0398854,-2.19178e-05,4.51192e-09,-3.28483e-13,-58451.7,-0.176787], Tmin=(666.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-480.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCCF) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFH) + group(CsCFFO) + polycyclic(s2_4_4_ane) + radical(cyclobutane) + radical(CsCsCsF1s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = 'F[C]=CC(F)=CO[C](F)F(13443)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {1,S} {7,D} {8,S}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u0 p0 c0 {6,S} {10,D} {12,S}
9  C u1 p0 c0 {2,S} {3,S} {5,S}
10 C u1 p0 c0 {4,S} {8,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-432.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2995,3025,975,1000,1300,1375,400,500,1630,1680,493,600,700,1144,1293,167,640,1190,246.828,247.534,248.015,248.482],'cm^-1')),
        HinderedRotor(inertia=(0.439683,'amu*angstrom^2'), symmetry=1, barrier=(18.9829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.435062,'amu*angstrom^2'), symmetry=1, barrier=(18.9863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434727,'amu*angstrom^2'), symmetry=1, barrier=(18.9699,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.045862,0.087391,-0.000105818,6.17521e-08,-1.40729e-11,-51928.8,30.4822], Tmin=(100,'K'), Tmax=(1073.43,'K')), NASAPolynomial(coeffs=[18.5411,0.0184692,-9.5052e-06,1.93429e-09,-1.41135e-13,-55899.4,-60.0567], Tmin=(1073.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-432.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFFHO) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(CdCFH) + radical(Csj(F1s)(F1s)(O2s-Cd)) + radical(Cdj(Cd-CdH)(F1s))"""),
)

species(
    label = 'F[C]1[CH]C(F)=C1(615)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 C u1 p0 c0 {4,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (86.8055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,346,659,817,1284,271,519,563,612,1379,180,180,1104.81,1106.27,1801.34],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06134,0.0383117,-3.38918e-05,1.45559e-08,-2.45109e-12,10513.6,15.9062], Tmin=(100,'K'), Tmax=(1431.45,'K')), NASAPolynomial(coeffs=[12.0907,0.0102863,-4.52462e-06,8.78888e-10,-6.24417e-14,7642.3,-36.0769], Tmin=(1431.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.8055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CdHH)(F1s)(Cd-CdH)_ring) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring)"""),
)

species(
    label = '[O][C](F)F(2580)',
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
    label = 'F[C](F)OC1=C(F)C=C1F(13444)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,D} {8,S}
7  C u0 p0 c0 {1,S} {6,D} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {11,S}
10 C u1 p0 c0 {3,S} {4,S} {5,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-369.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([217,343,449,587,660,812,790,914,798,948,2950,1000,493,600,700,1144,1293,189.999,191.218,193.893,197.185,198.068,1570.71,1572.43,1572.98],'cm^-1')),
        HinderedRotor(inertia=(0.00465433,'amu*angstrom^2'), symmetry=1, barrier=(8.13269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768517,'amu*angstrom^2'), symmetry=1, barrier=(21.5596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753638,0.0774401,-9.79806e-05,6.32529e-08,-1.64819e-11,-44369.6,26.3008], Tmin=(100,'K'), Tmax=(927.454,'K')), NASAPolynomial(coeffs=[12.7933,0.0255156,-1.40036e-05,2.89051e-09,-2.11282e-13,-46602.9,-30.8772], Tmin=(927.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFFHO) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F)) + radical(Csj(F1s)(F1s)(O2s-Cd)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C](F)OC1C(F)=C=C1F(13445)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {10,D}
8  C u0 p0 c0 {2,S} {6,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 C u0 p0 c0 {7,D} {8,D}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-266.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,104,186,245,407,330,466,761,907,1218,1388,493,600,700,1144,1293,217.002,217.019,217.076,1647.47,1647.49,1647.49,1647.5,1647.51],'cm^-1')),
        HinderedRotor(inertia=(0.19181,'amu*angstrom^2'), symmetry=1, barrier=(6.40915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191806,'amu*angstrom^2'), symmetry=1, barrier=(6.40962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50331,0.0882959,-0.000157938,1.51164e-07,-5.61753e-11,-31975.4,29.8158], Tmin=(100,'K'), Tmax=(805.412,'K')), NASAPolynomial(coeffs=[6.82931,0.0355314,-1.99125e-05,4.00836e-09,-2.83453e-13,-32302.1,4.96379], Tmin=(805.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(CsFFHO) + group(CdCddCF) + group(CdCddCF) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
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
    label = 'F[C](F)OC1=[C]C=C1F(13446)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {9,D}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {4,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (6.73652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,518,736,852,873,2950,1000,493,600,700,1144,1293,180,180,180,180,1099.13,1125.24,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.182275,'amu*angstrom^2'), symmetry=1, barrier=(4.19085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182275,'amu*angstrom^2'), symmetry=1, barrier=(4.19085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594263,0.0720641,-8.77302e-05,5.05609e-08,-1.12106e-11,935.372,27.2172], Tmin=(100,'K'), Tmax=(1112.58,'K')), NASAPolynomial(coeffs=[17.6603,0.0107063,-5.00489e-06,9.90399e-10,-7.17343e-14,-2862.03,-56.937], Tmin=(1112.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.73652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsFFHO) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cd-Cd-Cd-Cd(F)) + radical(Csj(F1s)(F1s)(O2s-Cd)) + radical(cyclobutadiene-C1)"""),
)

species(
    label = 'F[C](F)OC1[C]=C=C1F(13447)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u1 p0 c0 {2,S} {3,S} {4,S}
8  C u1 p0 c0 {5,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (161.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,145,326,398,834,1303,493,600,700,1144,1293,180,180,180,180,949.147,1345.51,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.171405,'amu*angstrom^2'), symmetry=1, barrier=(3.94094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171405,'amu*angstrom^2'), symmetry=1, barrier=(3.94094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.941811,0.0789287,-0.000146528,1.43967e-07,-5.41126e-11,19581.2,29.0544], Tmin=(100,'K'), Tmax=(819.054,'K')), NASAPolynomial(coeffs=[5.27716,0.0334622,-1.87699e-05,3.76502e-09,-2.65102e-13,19685.9,13.9787], Tmin=(819.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(CsFFHO) + group(CdCddCF) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13) + radical(Csj(F1s)(F1s)(O2s-Cs)) + radical(Cds_S)"""),
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
    label = 'F[C](F)OC1=C(F)[CH]C1F(13448)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {6,S} {9,D}
8  C u1 p0 c0 {6,S} {9,S} {12,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 C u1 p0 c0 {3,S} {4,S} {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-495.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,2950,1000,271,519,563,612,1379,493,600,700,1144,1293,180,180,180,180,1077.4,1077.54,1077.61,1078.16],'cm^-1')),
        HinderedRotor(inertia=(0.124935,'amu*angstrom^2'), symmetry=1, barrier=(2.87249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218741,'amu*angstrom^2'), symmetry=1, barrier=(18.0093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524719,0.074337,-7.64782e-05,3.80407e-08,-7.46528e-12,-59521.3,29.1479], Tmin=(100,'K'), Tmax=(1230.45,'K')), NASAPolynomial(coeffs=[17.3008,0.0198005,-9.99437e-06,2.01917e-09,-1.46486e-13,-63649.7,-55.266], Tmin=(1230.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-495.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(CsFFHO) + group(Cds-CdsCsOs) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring) + radical(Csj(F1s)(F1s)(O2s-Cd)) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C](F)OC1C(F)=[C]C1F(13449)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
8  C u0 p0 c0 {2,S} {6,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 C u1 p0 c0 {7,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-394.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,164,312,561,654,898,1207,1299,3167,246,474,533,1155,493,600,700,1144,1293,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01421,0.070178,-7.38543e-05,3.98793e-08,-8.83143e-12,-47360.6,28.7277], Tmin=(100,'K'), Tmax=(1071.27,'K')), NASAPolynomial(coeffs=[12.3626,0.0278045,-1.4523e-05,2.95678e-09,-2.14944e-13,-49792,-26.8033], Tmin=(1071.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CsFFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(Csj(F1s)(F1s)(O2s-Cs)) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
)

species(
    label = 'F[C]1[CH]C(F)=C1OC(F)F(13450)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
7  C u0 p0 c0 {5,S} {8,S} {9,D}
8  C u1 p0 c0 {3,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {7,D} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {12,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-572.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.233454,0.0798492,-8.7102e-05,4.58547e-08,-9.43841e-12,-68705.4,28.0986], Tmin=(100,'K'), Tmax=(1184.1,'K')), NASAPolynomial(coeffs=[18.5377,0.0180146,-8.76941e-06,1.75146e-09,-1.26689e-13,-73040.1,-63.3017], Tmin=(1184.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-572.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCCFH) + group(Cs-(Cds-Cds)CsHH) + group(CsFFHO) + group(Cds-CdsCsOs) + group(CdCsCdF) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Csj(Cs-F1sCdH)(Cd-CdF1s)(H)_ring) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
)

species(
    label = 'F[C]1[C]=C(F)C1OC(F)F(13451)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
8  C u1 p0 c0 {3,S} {6,S} {10,S}
9  C u0 p0 c0 {4,S} {6,S} {10,D}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-468.504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,519,593,830,1115,1166,1389,1413,3114,346,659,817,1284,246,474,533,1155,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634294,0.0707846,-6.85598e-05,3.21547e-08,-5.96753e-12,-56224.3,29.1016], Tmin=(100,'K'), Tmax=(1294.72,'K')), NASAPolynomial(coeffs=[16.9983,0.0202272,-9.98494e-06,1.99308e-09,-1.43404e-13,-60461.5,-54.0717], Tmin=(1294.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CsFFHO) + group(CdCsCdF) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(Cdj(Cs-CsF1sH)(Cd-CsF1s)_ring)"""),
)

species(
    label = 'F[C](F)OC1[C]=CC1(F)F(13452)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {7,S} {10,D} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {5,S}
10 C u1 p0 c0 {6,S} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-426.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,274,345,380,539,705,1166,1213,493,600,700,1144,1293,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995108,0.0711168,-7.7042e-05,4.32576e-08,-9.98035e-12,-51142.2,28.8096], Tmin=(100,'K'), Tmax=(1030.54,'K')), NASAPolynomial(coeffs=[12.004,0.0283864,-1.48461e-05,3.02255e-09,-2.19719e-13,-53411.2,-24.6331], Tmin=(1030.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-426.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFF) + group(CsFFHO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cd-Cs-Cs(F)(F)-Cd) + radical(Csj(F1s)(F1s)(O2s-Cs)) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'F[C]1C=[C]C1OC(F)(F)F(13453)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
8  C u1 p0 c0 {4,S} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {12,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-505.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,431,527,613,668,835,1199,1245,1322,346,659,817,1284,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666599,0.0670955,-6.10914e-05,2.67153e-08,-4.60027e-12,-60692.7,30.0913], Tmin=(100,'K'), Tmax=(1395.39,'K')), NASAPolynomial(coeffs=[17.6426,0.0184322,-8.77975e-06,1.72264e-09,-1.22529e-13,-65430.4,-57.464], Tmin=(1395.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-505.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCCFH) + group(CsFFFO) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cs(F)-Cs-Cd-Cd) + radical(CsCdCsF1s) + radical(cyclobutene-vinyl)"""),
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
    E0 = (-242.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (82.5744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (235.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-234.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-164.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-198.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-29.1452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-63.1048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-161.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-16.2042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-158.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (52.2473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (120.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (223.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (110.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (161.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (354.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (4.67177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-33.1639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (66.6233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-101.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-32.2137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (35.2588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-36.1407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['CF2O(49)', 'FC1=CC(F)=C1(307)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)O[CH]C1(F)C=C1F(13438)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]F(156)', '[O]C1[C](F)C=C1F(1459)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['FC1=CC2(F)C1OC2(F)F(13401)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['FC1=CC(F)=C1OC(F)F(13439)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['FC1=C=C(F)C1OC(F)F(13403)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(44.1897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['F[C](F)OC1C2(F)[CH]C12F(13440)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.44222e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c4_beta;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['F[C]1C2OC(F)(F)C1[C]2F(13441)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.96041e+10,'s^-1'), n=0.12, Ea=(179.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra;radadd_intra_cs] for rate rule [Rn2c4_beta;doublebond_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['F[C]1[CH]C2(F)C1OC2(F)F(13442)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.47579e+10,'s^-1'), n=0.224681, Ea=(80.8823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C]=CC(F)=CO[C](F)F(13443)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingle]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CF2O(49)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.33888e+13,'m^3/(mol*s)'), n=-2.06969, Ea=(94.2125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.21909771748945858, var=15.050572460742625, Tref=1000.0, N=2985, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C](F)F(2580)', 'FC1=CC(F)=C1(307)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'F[C](F)OC1=C(F)C=C1F(13444)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'F[C](F)OC1C(F)=C=C1F(13445)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.10362,'m^3/(mol*s)'), n=2.13572, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R_Ext-5R!H-R_Ext-2CS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Sp-2CS=1CCOSS_Ext-2CS-R_Ext-4R!H-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-4R!H-R_Ext-5R!H-R_Ext-2CS-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C](F)F(2580)', 'F[C]1[CH]C(F)=C1(615)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', 'F[C](F)OC1=[C]C=C1F(13446)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(157.445,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HF(38)', 'F[C](F)OC1[C]=C=C1F(13447)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(195.331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CF2(43)', '[O]C1[C](F)C=C1F(1459)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(6.60728,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C](F)OC1=C(F)[CH]C1F(13448)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.62e+13,'s^-1'), n=-0.14, Ea=(184.096,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_OneDe;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_OneDe/O;Cs_H_out_noH]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C](F)OC1C(F)=[C]C1F(13449)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S_cy4;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    products = ['F[C]1[CH]C(F)=C1OC(F)F(13450)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.95012e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C]1[C]=C(F)C1OC(F)F(13451)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.01114e+10,'s^-1'), n=0.8095, Ea=(157.607,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_noH]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C](F)OC1[C]=CC1(F)F(13452)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(182.661,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['F[C]1C=[C]C1OC(F)(F)F(13453)'],
    products = ['F[C](F)OC1[C](F)C=C1F(13400)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.94559e+07,'s^-1'), n=1.15307, Ea=(190.844,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF
Multiplied by reaction path degeneracy 3.0"""),
)

network(
    label = 'PDepNetwork #3580',
    isomers = [
        'F[C](F)OC1[C](F)C=C1F(13400)',
    ],
    reactants = [
        ('CF2O(49)', 'FC1=CC(F)=C1(307)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #3580',
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

