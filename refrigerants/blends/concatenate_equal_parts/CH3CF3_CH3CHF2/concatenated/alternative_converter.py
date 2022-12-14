import numpy as np
import pandas as pd
import os 
import re
import subprocess
import sys
import csv



file_name = 'chem.inp'

############### copies chem.inp file to fixed_chem.inp #############################


#copy chem.inp file into dups folder, and will then convert this copy into .cti
command = f'scp {file_name} fixed_{file_name}'
os.system(command)

############## read the lines of the fixed_chem.inp ############################


with open('fixed_chem.inp','r') as f:
    chemkin_file_lines = f.readlines()



############################# converts the dup_chem.inp files to .cti files #############################

counter = 0
while counter == 0: 
        
    #run the ck2cti command until the command passes or we get an error message that we this script doesn't account for
        ck2cti_command = subprocess.run(['ck2cti', f'--input=fixed_{file_name}',f'--transport=tran.dat'], capture_output=True)
        output = str(ck2cti_command.stdout)
        
        #if command passed without an error
        if re.search('PASSED',output):
            print('**************************command passed, converting to cti***************************')
            counter = 1 
            
        else: 
            
            #if the error has to do with undeclared duplicate reactions detected 
            if re.search('Undeclared duplicate reactions detected',output):
                
                #find the reaction numbers that the error message reports
                rxn_numbers = re.findall('Reaction ([0-9]+)', output)
                
                matched_lines = []
                
                #iterate through the chemkin file to find the lines where the matched reaction starts
                for rxn in rxn_numbers: 
                    for idx, chemkin_line in enumerate(chemkin_file_lines): 
            
                        if re.search(f'Chemkin #{rxn}', chemkin_line):
                            matched_lines.append(idx)
                        
                #now that we have the line numbers, let's find out where the next blank line is and add in a "DUPLICATE"
                for line in matched_lines: 
                    
                    iterations = 0
                    
                    count = 0 
                    while count == 0: 
                        
                        print(iterations)
                        start = line + iterations
                        print(chemkin_file_lines[start])

                        if re.search('^\n$',chemkin_file_lines[start]):
                            print('entered first loop')
                            chemkin_file_lines.insert(start, "DUPLICATE")
                            print(chemkin_file_lines[line+1])
                            count = 1

                        else: 
                            print('entered second loop')
                            iterations += 1
                            
                
                with open(f'fixed_{file_name}', 'w+') as f: 
                    for line in chemkin_file_lines: 
                        f.write(line)
                          
                    
            else:
                print('Another error occured!')
                print(output)
                