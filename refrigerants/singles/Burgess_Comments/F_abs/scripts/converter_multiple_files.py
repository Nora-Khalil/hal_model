"""
Converts Chemkin files with unmarked duplicates to Chemkin files with all duplicates marked.
Creates CTI file from "fixed" chemkin file. 
"""


import cantera as ct
import os 
import re
from subprocess import getoutput
import sys
import csv


print("Running Cantera Version: " + str(ct.__version__))


# Functions 

def copy_to_folder(file_name):
    '''copies Chemkin files to "copies" folder '''

    #copy folders so i dont screw up the original, and change into the new directory with copies
    os.makedirs('copies', exist_ok=True)

    #copy chem.inp file into dups folder, and will then convert this copy into .cti
    os.command = f'scp {file_name} copies/copy_{file_name}'
    os.system(os.command)
    os.system('scp tran.dat copies/tran.dat')


    
def convert_to_cti(file_name):
    '''Edits copy of Chemkin file so the duplicates are marked, then converts to CTI''' 

    x = 0
    while x == 0: 
        
        #"output" is the string of the output when I try to convert this file to a cti file. Will probably produce an error
        output = getoutput( f'ck2cti --input=copy_{file_name} --transport=tran.dat') 
        
        #if command passed without an error
        if re.search('PASSED',output):
            print('************************** command passed, converting to cti ***************************')
            x += 1 
            
        #if command generated the Duplicate error
        else:
            print('******************************* needs some work ****************************************')
            
            if re.search('Encountered\sunmarked\sduplicate\sreaction',output):
                
                match = re.search('See\slines\s([0-9]+)\sand\s([0-9]+)\sof\sthe\sinput\sfile',output)
                #capture the line numbers with the duplicate reactions
                
                line_numbers = [int(match.group(1)), int(match.group(2))]
                print(f'Unmarked duplicates on lines {match.group(1)} and {match.group(2)}')
                print('Editing chemkin file to allow conversion to .cti')
                
                #write the lines of the chemkin input file to a list so that I can insert the "DUPLICATE" statement
                with open(f'copy_{file_name}','r') as f:
                    data = f.readlines()
                    print(data[line_numbers[0]-1], data[line_numbers[1]-1])

                #start editing the .inp file below

                #'adjustments' will make sure that, even when I add an element in 'data', my index will still be correct
                adjustments = [0,1]
                for i,adjust in zip(line_numbers,adjustments): 
                    start = i+adjust-1
                    count = 0 
                    while count == 0: 
                        #if you don't see a blank line after the duplicated reaction line, keep going until you do
                        if not re.search('^\n', data[start]): 
                            print('no match')
                            print(start)
                            print(data[start])
                            start += 1
                        #when we get to the blank line after the reaction block, insert "DUPLICATE" and stop the loop for this line number
                        else: 
                            print('there is a match')
                            data.insert(start,'DUPLICATE')
                            count = 1 
                #now overwrite the input file with the change 
                with open(f'copy_{file_name}','w+') as f: 
                    for l in data: 
                        f.write(l)
                        x==0

            #if the command generated an error that is not the Duplicate error
            else:
                #if code ever gets to here, just cry
                print('There is another error, see Output')
                print(output)
                x += 1




##############################   execute code below.  ###################################################


#user defined variables
full_path = '/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/F_abs/DF_models/rmg_combustion_paper_FFCM/suppressants'
files_to_exclude = ['2-BTP', 'C3F7H', 'CF2BrCHCH2', 'CF2ClBr', 'CF3CBrCCF3', 'CF3CBrCF2', 'CF3CBrCHF', 'CF3CHCl2', 'CF3H', 'CHF2CBrCH2', 'halomethanes', 'HCFO-1233xf', 'HFO-1234yf',  'HFO-1336mzz']
unmarked_duplicates_file = 'chem_annotated.inp'

species_folders = [file for file in os.listdir(full_path) if file not in files_to_exclude]

for species in species_folders: 
    
    #switch to chemkin folder of each species
    chemkin_folder = os.path.join(full_path, species, 'chemkin')
    if os.path.exists(chemkin_folder):
        os.chdir(chemkin_folder)
    else: 
        continue 
        
    copy_to_folder(unmarked_duplicates_file)
    
    #switch to "copies" folder
    os.chdir('./copies')
    convert_to_cti(unmarked_duplicates_file)
    
    print(f'Completed {species}')
    

