import cantera as ct
import numpy as np
import pandas as pd
import os 
import re
from subprocess import getoutput
import sys
import csv


print("Running Cantera Version: " + str(ct.__version__))


'''copies chemkin files to dups folder '''


directory ='.'

#list all the input files in the chemkin folder
list_of_inp_files = [file for file in os.listdir(os.path.join(directory,'chemkin')) if re.match('chem0([0-9]+)\.inp', file) ]

#os.system('source activate ct_env') #if on local, this is cantera 2.6 beta
os.system('source activate cantera_env') #if on discovery (will be cantera 2.5)

#os.system('cd chemkin')
#os.getcwd()

os.chdir('./chemkin')

#copy folders so i dont screw up the originals, and change into the new directory with copies
os.makedirs('dups', exist_ok=True)

for file in list_of_inp_files: 
    os.command = f'scp {file} dups/dup_{file}'
    os.system(os.command)
    
dup_files = [file for file in os.listdir('./dups')]
os.system('scp tran.dat ./dups/tran.dat')
#os.system('cd dups') 

os.chdir('./dups')





'''converts all chemkin files to .cti files'''

#dup_files= ['dup_chem0156.inp'] #test with one file first 

for file in dup_files: 
    x = 0
    while x == 0: 
    #this is a string of the output when I try to convert this file to a cti file, will probably produce an error
        output = getoutput( f'ck2cti --input={file} --transport=tran.dat') 
        #if command passed without an error
        if re.search('PASSED',output):
            print('**************************command passed, converting to cti***************************')
            x += 1 
        #if command generated the Duplicate error
        else:
            print('*******************************needs some work****************************************')
            if re.search('Encountered\sunmarked\sduplicate\sreaction',output):
                match = re.search('See\slines\s([0-9]+)\sand\s([0-9]+)\sof\sthe\sinput\sfile',output)
                #capture the line numbers with the duplicate reactions
                line_numbers = [int(match.group(1)), int(match.group(2))]
                print(f'Unmarked duplicates on lines {match.group(1)} and {match.group(2)}')
                print('Editing chemkin file to allow conversion to .cti')
                #write the lines of the chemkin input file to a list so that I can insert the "DUPLICATE" statement
                with open(file,'r') as f:
                    data = f.readlines()
                    print(data[line_numbers[0]-1], data[line_numbers[1]-1])

                #start editing below

                #'adjustments' will make sure that, even when I add an element in 'data', my index will still be correct
                adjustments = [0,1]
                for i,adjust in zip(line_numbers,adjustments): 
                    start = i+adjust-1
                    count = 0 
                    while count == 0: 
                        #if you don't see a blank line after the duplicated reaction line, keep going until you do
                        if not re.search('^\\n', data[start]): 
                            print('no match')
                            print(start)
                            print(data[start])
                            start += 1
                        #when we get to the blank line after the reaction block, insert "DUPLICATE" and stop the loop for this line number
                        else: 
                            print('there is a match')
                            data.insert(start,'DUPLICATE')
                            count = 1 
                #now overwrite the file with the change 
                with open(f'{file}','w+') as f: 
                    for l in data: 
                        f.write(l)
                        x==0

            #if the command generated an error that is not the Duplicate error
            else:
                #if code ever gets to here, just cry
                print('There is another error, see Output')
                print(output)
                x += 1


'''Run this block if you want to see if all inp files were converted to cti files'''


cti_files = [file for file in os.listdir('.') if re.match('dup_chem0([0-9]+)\.cti', file)]

inp_files = [file for file in os.listdir('.') if re.match('dup_chem0([0-9]+)\.inp', file)]


print(f'Number of .cti files: {len(cti_files)} | Number of chemkin files: {len(inp_files)}')

if len(cti_files) == len(inp_files): 
    print('\nAll chemkin files converted')
else: 
    print('\nNot all chemkin files converted')



'''calculates flamespeeds for each .cti files at equivalence ratio = 1. Uses initial guess from the previous model'''

#os.chdir('/Users/nora/Code/projects/halogens/refrigerants/blends/C2H5F_CH2F2') #put slurm_task id as C2H5F_CH2F2 when I want to generalize the file?

cti_files = [file for file in os.listdir('.') if re.match('dup_chem0([0-9]+)\.cti', file)]

#cti_files = ['dup_chem0156.cti']
           
To = 298
Po = 1e5 # ct.one_atm
#vol_frac_list = np.arange(0.025, 0.25, step=0.01)
vol_frac_list =[.050]


header = ['species']
List_to_write_to_csv = []

os.chdir('../../')

for file in cti_files: 
    #make directory to store flame speed calculations, and a header for csv files
    d = f'./flame_calcs_different_models/{file}_CALC'
    os.makedirs(d, exist_ok=True)
    match1 = re.match('dup_chem0([0-9]+)\.cti',file) 
    gas = ct.Solution(f'./chemkin/dups/{file}')
    results = {}    
    header = ['species']
    #this part calculates the flame calcs for all of the vol_frac list, using the previous vol frac count as a guess for the current one
    for i in  range(len(vol_frac_list)):
        try:
           
            string = f'****************************starting new volume fraction: {vol_frac_list[i]}**************************'
            print(string)
           
            x = vol_frac_list[i]
            norm_ox = (1-x)*.21
            vol_frac_dict = {'CH2F2(1)': (x/3/norm_ox), 'C2H5F(2)': (x/3/norm_ox), 'C2H3F3(3)': (x/3/norm_ox), 'O2(4)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
            gas.TPX = To, Po, vol_frac_dict
            width = 0.08
            flame = ct.FreeFlame(gas, width=width)
            flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
            flame.max_time_step_count = 900
            loglevel = 1 
            if i!=0:
                d = f'./flame_calcs_different_models/{file}_CALC/test_{vol_frac_list[i-1]}.csv'
                if os.path.exists(d):  
                    arr2 = ct.SolutionArray(gas)
                    arr2.read_csv(d)
                    flame.set_initial_guess(data=arr2)
                    print(' initial guess has been set')
            #"False" stops the calculation from retrying over and over, thanks Chao 
            flame.solve(loglevel=loglevel, auto=False)
            #flame.solve(loglevel=loglevel, auto=True)
            Su = flame.velocity[0]
            results[x] = Su
            sltn = flame.to_solution_array()
            df1 = sltn.to_pandas()
            df1.to_csv(f'./flame_calcs_different_models/{file}_CALC/test_{x}.csv')
            
            

             #delete first column of csv file to avoid error
            df2 = pd.read_csv(f'./flame_calcs_different_models/{file}_CALC/test_{x}.csv')
            first_column = df2.columns[0]
            # Delete first
            df2 = df2.drop([first_column], axis=1)
            df2.to_csv(f'./flame_calcs_different_models/{file}_CALC/test_{x}.csv', index=False)
        except Exception as e: 
            print(f'********************passed volume fraction:{vol_frac_list[i]}, error: {e}*************************************')
            pass

    vol_fracs = list(results.keys())
    flame_speeds = list(results.values())


    print("volume fractions are:")
    print(vol_fracs)
 
    print("flame speeds are:")
    print(flame_speeds)
    
    header.append(vol_fracs)
    header.append([f'speed_for_volfrac{x}' for x in flame_speeds])
    
    values = [match1.group(1),vol_fracs, flame_speeds]
    
    List_to_write_to_csv.append(values)
    
with open('final_calcs.csv', 'w+') as g:
    writers = csv.writer(g)
    writers.writerow(header)
    for i in List_to_write_to_csv:
    	writers.writerow(i)
