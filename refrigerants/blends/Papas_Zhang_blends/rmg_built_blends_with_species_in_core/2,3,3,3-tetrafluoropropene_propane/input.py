 

thermolibs = [
'primaryThermoLibrary',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'Fluorine',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes'
]

thermolibs_Creg = [
'primaryThermoLibrary',
'FFCM1(-)',
'DFT_QCI_thermo',
]



database(
thermoLibraries = thermolibs,
reactionLibraries = ['FFCM1(-)','halogens_pdep'],
seedMechanisms = ['FFCM1(-)'],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)


    
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)
    
    
    


species(
    label = '2,3,3,3-tetrafluoropropene',
    reactive = True,
    structure = SMILES('C=C(F)C(F)(F)F')
)      
        
        

species(
    label = 'OH',
    reactive = True,
    structure = SMILES('[OH]')
)      
        
        

species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)      
        
        

species(
    label = 'H2O',
    reactive = True,
    structure = SMILES('O')
)      
        
        

species(
    label = 'H',
    reactive = True,
    structure = SMILES('[H]')
)      
        
        

species(
    label = 'O',
    reactive = True,
    structure = SMILES('[O]')
)      
        
        

species(
    label = 'H2',
    reactive = True,
    structure = SMILES('[H][H]')
)      
        
        

species(
    label = 'Ar',
    reactive = True,
    structure = SMILES('[Ar]')
)      
        
        

species(
    label = 'HO2',
    reactive = True,
    structure = SMILES('[O]O')
)      
        
        

species(
    label = 'He',
    reactive = True,
    structure = SMILES('[He]')
)      
        
        

species(
    label = 'H2O2',
    reactive = True,
    structure = SMILES('OO')
)      
        
        

species(
    label = 'CO',
    reactive = True,
    structure = SMILES('[C-]#[O+]')
)      
        
        

species(
    label = 'CO2',
    reactive = True,
    structure = SMILES('O=C=O')
)      
        
        

species(
    label = 'HCO',
    reactive = True,
    structure = SMILES('[CH]=O')
)      
        
        

species(
    label = 'C(T)',
    reactive = True,
    structure = SMILES('[C]')
)      
        
        

species(
    label = 'CH',
    reactive = True,
    structure = SMILES('[CH]')
)      
        
        

species(
    label = 'CH2(T)',
    reactive = True,
    structure = SMILES('[CH2]')
)      
        
        

species(
    label = 'CH3',
    reactive = True,
    structure = SMILES('[CH3]')
)      
        
        

species(
    label = 'CH2O',
    reactive = True,
    structure = SMILES('C=O')
)      
        
        

species(
    label = 'HCCO',
    reactive = True,
    structure = SMILES('C#C[O]')
)      
        
        

species(
    label = 'C2H',
    reactive = True,
    structure = SMILES('[C]#C')
)      
        
        

species(
    label = 'C2H2',
    reactive = True,
    structure = SMILES('C#C')
)      
        
        

species(
    label = 'H2CC',
    reactive = True,
    structure = SMILES('[C]=C')
)      
        
        

species(
    label = 'CH3OH',
    reactive = True,
    structure = SMILES('CO')
)      
        
        

species(
    label = 'CH3O',
    reactive = True,
    structure = SMILES('C[O]')
)      
        
        

species(
    label = 'CH2CO',
    reactive = True,
    structure = SMILES('C=C=O')
)      
        
        

species(
    label = 'C2H3',
    reactive = True,
    structure = SMILES('[CH]=C')
)      
        
        

species(
    label = 'C2H4',
    reactive = True,
    structure = SMILES('C=C')
)      
        
        

species(
    label = 'CH4',
    reactive = True,
    structure = SMILES('C')
)      
        
        

species(
    label = 'C2H6',
    reactive = True,
    structure = SMILES('CC')
)      
        
        

species(
    label = 'C2H5',
    reactive = True,
    structure = SMILES('C[CH2]')
)      
        
        

species(
    label = 'CH2OH',
    reactive = True,
    structure = SMILES('[CH2]O')
)      
        
        

species(
    label = 'CH3CO',
    reactive = True,
    structure = SMILES('C[C]=O')
)      
        
        

species(
    label = 'CH2CHO',
    reactive = True,
    structure = SMILES('[CH2]C=O')
)      
        
        

species(
    label = 'CH3CHO',
    reactive = True,
    structure = SMILES('CC=O')
)      
        
        

species(
    label = 'F',
    reactive = True,
    structure = SMILES('[F]')
)      
        
        

species(
    label = 'HF',
    reactive = True,
    structure = SMILES('F')
)      
        
        

species(
    label = 'CH2F2',
    reactive = True,
    structure = SMILES('FCF')
)      
        
        

species(
    label = 'CHF3',
    reactive = True,
    structure = SMILES('FC(F)F')
)      
        
        

species(
    label = 'CF2',
    reactive = True,
    structure = SMILES('F[C]F')
)      
        
        

species(
    label = 'CF4',
    reactive = True,
    structure = SMILES('FC(F)(F)F')
)      
        
        

species(
    label = 'CF3',
    reactive = True,
    structure = SMILES('F[C](F)F')
)      
        
        

species(
    label = 'CH2F',
    reactive = True,
    structure = SMILES('[CH2]F')
)      
        
        

species(
    label = 'CHFO',
    reactive = True,
    structure = SMILES('O=CF')
)      
        
        

species(
    label = 'CF3O',
    reactive = True,
    structure = SMILES('[O]C(F)(F)F')
)      
        
        

species(
    label = 'CF2O',
    reactive = True,
    structure = SMILES('O=C(F)F')
)      
        
        

species(
    label = 'CFO',
    reactive = True,
    structure = SMILES('O=[C]F')
)      
        
        

species(
    label = 'CF3-CF3',
    reactive = True,
    structure = SMILES('FC(F)(F)C(F)(F)F')
)      
        
        

species(
    label = 'CF3-CHF',
    reactive = True,
    structure = SMILES('F[CH]C(F)(F)F')
)      
        
        

species(
    label = 'CHFCF2',
    reactive = True,
    structure = SMILES('FC=C(F)F')
)      
        
        

species(
    label = 'CH2CHF',
    reactive = True,
    structure = SMILES('C=CF')
)      
        
        

species(
    label = 'CH2CF2',
    reactive = True,
    structure = SMILES('C=C(F)F')
)      
        
        

species(
    label = 'C2HF',
    reactive = True,
    structure = SMILES('C#CF')
)      
        
        

species(
    label = 'CF2CF2',
    reactive = True,
    structure = SMILES('FC(F)=C(F)F')
)      
        
        

species(
    label = 'CHF2-CHF',
    reactive = True,
    structure = SMILES('F[CH]C(F)F')
)      
        
        

species(
    label = 'CH2F-CF2',
    reactive = True,
    structure = SMILES('FC[C](F)F')
)      
        
        

species(
    label = 'CHFCF[Z]',
    reactive = True,
    structure = SMILES('F[C]=CF')
)      
        
        

species(
    label = 'CF2CH',
    reactive = True,
    structure = SMILES('[CH]=C(F)F')
)      
        
        

species(
    label = 'CF3COF',
    reactive = True,
    structure = SMILES('O=C(F)C(F)(F)F')
)      
        
        

species(
    label = 'CH2CFO',
    reactive = True,
    structure = SMILES('[CH2]C(=O)F')
)      
        
        

species(
    label = 'CHF2',
    reactive = True,
    structure = SMILES('F[CH]F')
)      
        
        

species(
    label = 'CF3CCH',
    reactive = True,
    structure = SMILES('C#CC(F)(F)F')
)      
        
        

species(
    label = 'E-CHCFCF3',
    reactive = True,
    structure = SMILES('[CH]=C(F)C(F)(F)F')
)      
        
        

species(
    label = 'CHOCFCF3',
    reactive = True,
    structure = SMILES('O=C[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'CH3CFCF3',
    reactive = True,
    structure = SMILES('C[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'CF3CFCO',
    reactive = True,
    structure = SMILES('O=C=C(F)C(F)(F)F')
)      
        
        

species(
    label = 'FCC(F)=C(F)F',
    reactive = True,
    structure = SMILES('FCC(F)=C(F)F')
)      
        
        

species(
    label = 'C=C(F)[C](F)F',
    reactive = True,
    structure = SMILES('C=C(F)[C](F)F')
)      
        
        

species(
    label = 'FC[C]C(F)(F)F',
    reactive = True,
    structure = SMILES('FC[C]C(F)(F)F')
)      
        
        

species(
    label = '[CH]C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH]C(F)C(F)(F)F')
)      
        
        

species(
    label = 'C=C=C(F)F',
    reactive = True,
    structure = SMILES('C=C=C(F)F')
)      
        
        

species(
    label = 'F[C]C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C]C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(C[C](F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(C[C](F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(O[O])C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(O[O])C(F)(F)F')
)      
        
        

species(
    label = 'FC=CC(F)(F)F',
    reactive = True,
    structure = SMILES('FC=CC(F)(F)F')
)      
        
        

species(
    label = 'FC(F)=CC(F)F',
    reactive = True,
    structure = SMILES('FC(F)=CC(F)F')
)      
        
        

species(
    label = 'F[C]CC(F)(F)F',
    reactive = True,
    structure = SMILES('F[C]CC(F)(F)F')
)      
        
        

species(
    label = 'FC=C=C(F)F',
    reactive = True,
    structure = SMILES('FC=C=C(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(F)F')
)      
        
        

species(
    label = 'O=O',
    reactive = True,
    structure = SMILES('O=O')
)      
        
        

species(
    label = '[CH2]C([O])(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C([O])(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)COO1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)COO1')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)CO1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)CO1')
)      
        
        

species(
    label = 'F[C](F)C[C](F)F',
    reactive = True,
    structure = SMILES('F[C](F)C[C](F)F')
)      
        
        

species(
    label = 'FC(F)[C]C(F)F',
    reactive = True,
    structure = SMILES('FC(F)[C]C(F)F')
)      
        
        

species(
    label = 'F[C](C[C](F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C[C](F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(C=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(C=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(=CC(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC(=CC(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]O[C](F)F',
    reactive = True,
    structure = SMILES('[O]O[C](F)F')
)      
        
        

species(
    label = '[O]C[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]C[C](F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OCC(F)=C(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)=C(F)F')
)      
        
        

species(
    label = 'F[C](F)C[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](F)C[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'F[CH]C(F)[C](F)F',
    reactive = True,
    structure = SMILES('F[CH]C(F)[C](F)F')
)      
        
        

species(
    label = 'C1OO1',
    reactive = True,
    structure = SMILES('C1OO1')
)      
        
        

species(
    label = 'CC(F)(O[O])C(F)(F)F',
    reactive = True,
    structure = SMILES('CC(F)(O[O])C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(OO)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(OO)C(F)(F)F')
)      
        
        

species(
    label = 'FC1(F)CO1',
    reactive = True,
    structure = SMILES('FC1(F)CO1')
)      
        
        

species(
    label = 'OC[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('OC[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'FC1(F)CC1(F)F',
    reactive = True,
    structure = SMILES('FC1(F)CC1(F)F')
)      
        
        

species(
    label = 'FC=CC(F)(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC=CC(F)(F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)([C](F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)([C](F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]C[O]',
    reactive = True,
    structure = SMILES('[O]C[O]')
)      
        
        

species(
    label = '[O]C=O',
    reactive = True,
    structure = SMILES('[O]C=O')
)      
        
        

species(
    label = 'O=C[C](F)F',
    reactive = True,
    structure = SMILES('O=C[C](F)F')
)      
        
        

species(
    label = 'FC=C(F)C(F)F',
    reactive = True,
    structure = SMILES('FC=C(F)C(F)F')
)      
        
        

species(
    label = 'F[CH][C](F)F',
    reactive = True,
    structure = SMILES('F[CH][C](F)F')
)      
        
        

species(
    label = 'FC1C(F)C1(F)F',
    reactive = True,
    structure = SMILES('FC1C(F)C1(F)F')
)      
        
        

species(
    label = 'CC(F)=C(F)F',
    reactive = True,
    structure = SMILES('CC(F)=C(F)F')
)      
        
        

species(
    label = '[O]OCC(F)(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)(F)C(F)(F)F')
)      
        
        

species(
    label = 'C=C(F)C(F)(F)O[O]',
    reactive = True,
    structure = SMILES('C=C(F)C(F)(F)O[O]')
)      
        
        

species(
    label = 'F[C]1COOC1(F)F',
    reactive = True,
    structure = SMILES('F[C]1COOC1(F)F')
)      
        
        

species(
    label = '[O]OCF',
    reactive = True,
    structure = SMILES('[O]OCF')
)      
        
        

species(
    label = 'F[C](CC(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](CC(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)F')
)      
        
        

species(
    label = 'O=CC(F)=C(F)F',
    reactive = True,
    structure = SMILES('O=CC(F)=C(F)F')
)      
        
        

species(
    label = 'FC([CH]C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC([CH]C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(OOC(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(OOC(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[C]=CC(F)(F)F',
    reactive = True,
    structure = SMILES('[C]=CC(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)C(F)(F)F')
)      
        
        

species(
    label = '[O]CC(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]CC(F)C(F)(F)F')
)      
        
        

species(
    label = 'O=CC(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=CC(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)=CC(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC(F)=CC(F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(CC(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(CC(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)[CH]O1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)[CH]O1')
)      
        
        

species(
    label = '[O]OCC(=O)F',
    reactive = True,
    structure = SMILES('[O]OCC(=O)F')
)      
        
        

species(
    label = 'O[CH]C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O[CH]C(F)C(F)(F)F')
)      
        
        

species(
    label = 'F[C](CC=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](CC=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[C]#CF',
    reactive = True,
    structure = SMILES('[C]#CF')
)      
        
        

species(
    label = '[C]=CF',
    reactive = True,
    structure = SMILES('[C]=CF')
)      
        
        

species(
    label = 'FC1=CC1(F)F',
    reactive = True,
    structure = SMILES('FC1=CC1(F)F')
)      
        
        

species(
    label = 'O=C=[C]F',
    reactive = True,
    structure = SMILES('O=C=[C]F')
)      
        
        

species(
    label = 'FC(F)[C]C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC(F)[C]C(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC1(F)OO1',
    reactive = True,
    structure = SMILES('FC1(F)OO1')
)      
        
        

species(
    label = 'F[C](C1CC1(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C1CC1(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OCC(F)(C=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)(C=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]C([O])(F)F',
    reactive = True,
    structure = SMILES('[O]C([O])(F)F')
)      
        
        

species(
    label = '[O]OC(F)(CC=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(CC=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[CH]=C[O]',
    reactive = True,
    structure = SMILES('[CH]=C[O]')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)[CH]CC(F)(C(F)(F)F)OO1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)[CH]CC(F)(C(F)(F)F)OO1')
)      
        
        

species(
    label = 'F[C](C1CC(F)(C(F)(F)F)OO1)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C1CC(F)(C(F)(F)F)OO1)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(C1CC(F)(C(F)(F)F)OO1)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(C1CC(F)(C(F)(F)F)OO1)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(C=O)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(C=O)C(F)(F)F')
)      
        
        

species(
    label = '[O]C1OOC1(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]C1OOC1(F)C(F)(F)F')
)      
        
        

species(
    label = 'OC=CF',
    reactive = True,
    structure = SMILES('OC=CF')
)      
        
        

species(
    label = '[O]OC(O)C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(O)C(F)C(F)(F)F')
)      
        
        

species(
    label = 'O=CCF',
    reactive = True,
    structure = SMILES('O=CCF')
)      
        
        

species(
    label = 'O=COO[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=COO[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'O=[C]C(F)[C](F)F',
    reactive = True,
    structure = SMILES('O=[C]C(F)[C](F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)CC1(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)CC1(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC=COC(F)(F)F',
    reactive = True,
    structure = SMILES('FC=COC(F)(F)F')
)      
        
        

species(
    label = 'FC[C]C(F)(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC[C]C(F)(F)C(F)(F)F')
)      
        
        

species(
    label = 'OO[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('OO[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)[CH]C(F)(C(F)(F)F)OOC1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)[CH]C(F)(C(F)(F)F)OOC1')
)      
        
        

species(
    label = 'F[C](C1OOCC1(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C1OOCC1(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'OOC(F)([C]1CC(F)(C(F)(F)F)OO1)C(F)(F)F',
    reactive = True,
    structure = SMILES('OOC(F)([C]1CC(F)(C(F)(F)F)OO1)C(F)(F)F')
)      
        
        

species(
    label = '[O]C(F)(CC(=O)C(F)(OO)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]C(F)(CC(=O)C(F)(OO)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'C=C([O])C(F)(OO)C(F)(F)F',
    reactive = True,
    structure = SMILES('C=C([O])C(F)(OO)C(F)(F)F')
)      
        
        

species(
    label = 'O=C(F)C(=O)CC(F)(F)F',
    reactive = True,
    structure = SMILES('O=C(F)C(=O)CC(F)(F)F')
)      
        
        

species(
    label = 'CC(=O)[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('CC(=O)[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'CC(=O)C(F)(O[O])C(F)(F)F',
    reactive = True,
    structure = SMILES('CC(=O)C(F)(O[O])C(F)(F)F')
)      
        
        

species(
    label = 'O=C1C(F)C1(F)F',
    reactive = True,
    structure = SMILES('O=C1C(F)C1(F)F')
)      
        
        

species(
    label = 'O=[C]C(=O)F',
    reactive = True,
    structure = SMILES('O=[C]C(=O)F')
)      
        
        

species(
    label = 'O=[C]C(F)(F)[CH]F',
    reactive = True,
    structure = SMILES('O=[C]C(F)(F)[CH]F')
)      
        
        

species(
    label = '[O]C1(F)OO1',
    reactive = True,
    structure = SMILES('[O]C1(F)OO1')
)      
        
        

species(
    label = 'propane',
    reactive = True,
    structure = SMILES('CCC')
)      
        
        

species(
    label = 'C[CH]C',
    reactive = True,
    structure = SMILES('C[CH]C')
)      
        
        

species(
    label = '[CH2]CC',
    reactive = True,
    structure = SMILES('[CH2]CC')
)      
        
        

species(
    label = 'CO[O]',
    reactive = True,
    structure = SMILES('CO[O]')
)      
        
        

species(
    label = 'CCO[O]',
    reactive = True,
    structure = SMILES('CCO[O]')
)      
        
        

species(
    label = 'CC[O]',
    reactive = True,
    structure = SMILES('CC[O]')
)      
        
        

species(
    label = '[CH]C-2',
    reactive = True,
    structure = SMILES('[CH]C')
)      
        
        

species(
    label = 'COO',
    reactive = True,
    structure = SMILES('COO')
)      
        
        

species(
    label = '[CH2]C[CH2]',
    reactive = True,
    structure = SMILES('[CH2]C[CH2]')
)      
        
        

species(
    label = 'C=CCC',
    reactive = True,
    structure = SMILES('C=CCC')
)      
        
        

species(
    label = '[CH2]CC=C',
    reactive = True,
    structure = SMILES('[CH2]CC=C')
)      
        
        

species(
    label = '[CH2]CO',
    reactive = True,
    structure = SMILES('[CH2]CO')
)      
        
        

species(
    label = 'C=CC',
    reactive = True,
    structure = SMILES('C=CC')
)      
        
        

species(
    label = 'CC(C)O[O]',
    reactive = True,
    structure = SMILES('CC(C)O[O]')
)      
        
        

species(
    label = 'CCCO[O]',
    reactive = True,
    structure = SMILES('CCCO[O]')
)      
        
        

species(
    label = '[CH2]CCOO',
    reactive = True,
    structure = SMILES('[CH2]CCOO')
)      
        
        

species(
    label = '[CH2]C=C',
    reactive = True,
    structure = SMILES('[CH2]C=C')
)      
        
        

species(
    label = 'C=[C]C',
    reactive = True,
    structure = SMILES('C=[C]C')
)      
        
        

species(
    label = '[CH]=CC',
    reactive = True,
    structure = SMILES('[CH]=CC')
)      
        
        

species(
    label = 'C[C]C-2',
    reactive = True,
    structure = SMILES('C[C]C')
)      
        
        

species(
    label = 'C[CH]CC',
    reactive = True,
    structure = SMILES('C[CH]CC')
)      
        
        

species(
    label = 'C=C[C]=O',
    reactive = True,
    structure = SMILES('C=C[C]=O')
)      
        
        

species(
    label = 'C=C=C',
    reactive = True,
    structure = SMILES('C=C=C')
)      
        
        

species(
    label = 'C=CCO[O]',
    reactive = True,
    structure = SMILES('C=CCO[O]')
)      
        
        

species(
    label = '[CH]1COOC1',
    reactive = True,
    structure = SMILES('[CH]1COOC1')
)      
        
        

species(
    label = 'C#CC=O',
    reactive = True,
    structure = SMILES('C#CC=O')
)      
        
        

species(
    label = 'C1=CCC1',
    reactive = True,
    structure = SMILES('C1=CCC1')
)      
        
        

species(
    label = 'C=CC=C',
    reactive = True,
    structure = SMILES('C=CC=C')
)      
        
        

species(
    label = 'C#CC[CH2]',
    reactive = True,
    structure = SMILES('C#CC[CH2]')
)      
        
        

species(
    label = '[CH2]C=C[CH2]',
    reactive = True,
    structure = SMILES('[CH2]C=C[CH2]')
)      
        
        

species(
    label = 'C1=CO1',
    reactive = True,
    structure = SMILES('C1=CO1')
)      
        
        

species(
    label = 'C=C[CH]C',
    reactive = True,
    structure = SMILES('C=C[CH]C')
)      
        
        

species(
    label = '[CH]=CCC',
    reactive = True,
    structure = SMILES('[CH]=CCC')
)      
        
        

species(
    label = 'C[C]CC',
    reactive = True,
    structure = SMILES('C[C]CC')
)      
        
        

species(
    label = 'C1=COOC1',
    reactive = True,
    structure = SMILES('C1=COOC1')
)      
        
        

species(
    label = '[CH]=C=C',
    reactive = True,
    structure = SMILES('[CH]=C=C')
)      
        
        

species(
    label = 'C=C=CO[O]',
    reactive = True,
    structure = SMILES('C=C=CO[O]')
)      
        
        

species(
    label = 'C#CCO[O]',
    reactive = True,
    structure = SMILES('C#CCO[O]')
)      
        
        

species(
    label = 'C=[C]C=O',
    reactive = True,
    structure = SMILES('C=[C]C=O')
)      
        
        

species(
    label = 'C#CC',
    reactive = True,
    structure = SMILES('C#CC')
)      
        
        

species(
    label = 'C=C=CC',
    reactive = True,
    structure = SMILES('C=C=CC')
)      
        
        

species(
    label = '[C]1=COOC1',
    reactive = True,
    structure = SMILES('[C]1=COOC1')
)      
        
        

species(
    label = 'C=[C]C1OO1',
    reactive = True,
    structure = SMILES('C=[C]C1OO1')
)      
        
        

species(
    label = 'C#C[CH]C',
    reactive = True,
    structure = SMILES('C#C[CH]C')
)      
        
        

species(
    label = 'C=COC=O',
    reactive = True,
    structure = SMILES('C=COC=O')
)      
        
        

species(
    label = '[C]=CC',
    reactive = True,
    structure = SMILES('[C]=CC')
)      
        
        

species(
    label = '[O]C=CC[O]',
    reactive = True,
    structure = SMILES('[O]C=CC[O]')
)      
        
        

species(
    label = '[O]C=CC=O',
    reactive = True,
    structure = SMILES('[O]C=CC=O')
)      
        
        

species(
    label = 'O=CCC=O',
    reactive = True,
    structure = SMILES('O=CCC=O')
)      
        
        

species(
    label = 'O=[C]CC=O',
    reactive = True,
    structure = SMILES('O=[C]CC=O')
)      
        
        

species(
    label = '[CH2]C(=O)C=O',
    reactive = True,
    structure = SMILES('[CH2]C(=O)C=O')
)      
        
        

species(
    label = 'C#CC(=C)[O]',
    reactive = True,
    structure = SMILES('C#CC(=C)[O]')
)      
        
        

species(
    label = 'C#C[C]=O',
    reactive = True,
    structure = SMILES('C#C[C]=O')
)      
        
        

species(
    label = 'C#CC=C',
    reactive = True,
    structure = SMILES('C#CC=C')
)      
        
        

species(
    label = 'C=[C]C=C',
    reactive = True,
    structure = SMILES('C=[C]C=C')
)      
        
        

species(
    label = 'C=C1[CH]O1',
    reactive = True,
    structure = SMILES('C=C1[CH]O1')
)      
        
        

species(
    label = 'C=C(C=O)O[O]',
    reactive = True,
    structure = SMILES('C=C(C=O)O[O]')
)      
        
        

species(
    label = '[CH2]C1(C=O)OO1',
    reactive = True,
    structure = SMILES('[CH2]C1(C=O)OO1')
)      
        
        

species(
    label = '[C]=CC=O',
    reactive = True,
    structure = SMILES('[C]=CC=O')
)      
        
        

species(
    label = 'C=C[C]C',
    reactive = True,
    structure = SMILES('C=C[C]C')
)      
        
        

species(
    label = '[CH2]C1CC1',
    reactive = True,
    structure = SMILES('[CH2]C1CC1')
)      
        
        

species(
    label = 'CC=CC',
    reactive = True,
    structure = SMILES('CC=CC')
)      
        
        

species(
    label = 'CC=CCO[O]',
    reactive = True,
    structure = SMILES('CC=CCO[O]')
)      
        
        

species(
    label = 'C=CC(C)O[O]',
    reactive = True,
    structure = SMILES('C=CC(C)O[O]')
)      
        
        

species(
    label = '[C]1=CCC1',
    reactive = True,
    structure = SMILES('[C]1=CCC1')
)      
        
        

species(
    label = '[C]1CCC1',
    reactive = True,
    structure = SMILES('[C]1CCC1')
)      
        
        

species(
    label = 'C=C=CCO[O]',
    reactive = True,
    structure = SMILES('C=C=CCO[O]')
)      
        
        

species(
    label = 'C#CC(C)O[O]',
    reactive = True,
    structure = SMILES('C#CC(C)O[O]')
)      
        
        

species(
    label = 'C#C[C]C',
    reactive = True,
    structure = SMILES('C#C[C]C')
)      
        
        

species(
    label = '[C]=CC=C-2',
    reactive = True,
    structure = SMILES('[C]=CC=C')
)      
        
        

species(
    label = '[CH]=C=CCOO',
    reactive = True,
    structure = SMILES('[CH]=C=CCOO')
)      
        
        

species(
    label = '[CH2]C=C=C=O',
    reactive = True,
    structure = SMILES('[CH2]C=C=C=O')
)      
        
        

species(
    label = 'C#CC1CO1',
    reactive = True,
    structure = SMILES('C#CC1CO1')
)      
        
        

species(
    label = 'CC1[C]=COO1',
    reactive = True,
    structure = SMILES('CC1[C]=COO1')
)      
        
        

species(
    label = 'C=C1[CH]C1',
    reactive = True,
    structure = SMILES('C=C1[CH]C1')
)      
        
        

species(
    label = '[CH]=C1CC1',
    reactive = True,
    structure = SMILES('[CH]=C1CC1')
)      
        
        

species(
    label = '[CH]=C=C1CO1',
    reactive = True,
    structure = SMILES('[CH]=C=C1CO1')
)      
        
        

species(
    label = 'O=C=C1[CH]C1',
    reactive = True,
    structure = SMILES('O=C=C1[CH]C1')
)      
        
        

species(
    label = 'C=CC(=O)[C]=O',
    reactive = True,
    structure = SMILES('C=CC(=O)[C]=O')
)      
        
        
    
simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)


generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=8,
    maximumOxygenAtoms=6,
    #maximumHeavyAtoms=24,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = True,
)

options(
    units = "si",
    generateSeedEachIteration = True,
    generateOutputHTML = True,
    generatePlots = True,
    saveSimulationProfiles = True,
    saveEdgeSpecies = False,
    keepIrreversible = True,
    verboseComments = False,
)
    
    
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure= [(1.0,'bar'),(10.0,'bar')],
        nSims=10,
        initialMoleFractions={
        "2,3,3,3-tetrafluoropropene": 0.5,
        "propane": 0.5,
        "O2": 1,
        "N2": 3.76,

        },
        # terminationConversion={
        # 'halogen': 0.999,
        # },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        #sensitivity=['halogen','OH'],
        #sensitivityThreshold=0.001,
        )
        
model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

    
    
    