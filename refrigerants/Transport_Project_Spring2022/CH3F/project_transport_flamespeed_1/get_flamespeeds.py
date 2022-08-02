import os
import cantera as ct
import re 

print("Running Cantera Version: " + str(ct.__version__))

flamespeeds = []
species = []

list_of_cti_files = os.listdir('../project_transport_cantera_1')
print(list_of_cti_files)

for cti in list_of_cti_files: 
    match = re.match('chem([0-9]+)\.cti', cti)
    gas = ct.Solution(f'../project_transport_cantera_1/chem{match.group(1)}.cti')

    To = 298
    Po = 1e5 # ct.one_atm

    #only testing one equivalence ratio = 1 
    volume_fraction_list = [0.123]

    results = {}

    for x in volume_fraction_list: 
        norm_ox = (1-x)*.21
        volume_fraction_dict = {'CH3F(1)': (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
        gas.TPX = To, Po, volume_fraction_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        flame.solve(loglevel=loglevel, auto=True)
        Su = flame.u[0]
        results[x] = Su
        flamespeeds.append(Su)
        species.append(match.group(1))
        #save to csv file
        directory = './data'
        with open(os.path.join(directory,f'flamespeeds{match.group(1)}.csv'),'w') as f:
            f.write("flamespeeds:")
            f.write(str(flamespeeds))
            f.write("species:") 
            f.write(str(species))

#save to csv file
directory = './data'
with open(os.path.join(directory,'final.csv'),'w') as f:
    f.write("flamespeeds:")
    f.write(str(flamespeeds))
    f.write("species:") 
    f.write(str(species))







