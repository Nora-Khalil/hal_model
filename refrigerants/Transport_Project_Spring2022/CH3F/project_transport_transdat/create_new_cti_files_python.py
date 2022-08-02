import os
import cantera as ct 
import re


list_of_directories = os.listdir('.')
list_of_trandat_files = [file for file in list_of_directories if re.match('tran[0-9]+\.dat', file)] #list of only new .dat files 

os.system('source activate cantera_env')

for file in list_of_trandat_files: 
    match = re.match('tran([0-9]+)\.dat', file)
    os_command1 = f'ck2cti --input=chem.inp --transport={file}'
    os.system(os_command1)
    os_command2 = f'mv chem.cti chem{match.group(1)}\.cti'
    os.system(os_command2)
    
    
list_of_directories = os.listdir('.')
list_of_cti_files = [file for file in list_of_directories if re.match('chem[0-9]+\.cti', file)]

#move to a new folder
for cti in list_of_cti_files:
       match = re.match('chem([0-9]+)\.cti',cti)
       os_command3 = f'mv chem{match.group(1)}.cti ../project_transport_cantera'
       os.system(os_command3)

