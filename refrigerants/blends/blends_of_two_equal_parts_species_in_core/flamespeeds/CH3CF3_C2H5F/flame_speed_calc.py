

from os import uname
import cantera as ct
import numpy as np
import pandas as pd
import os
import csv

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

    
directory = '/work/westgroup/nora/Code/projects/halogens/refrigerants/blends/blends_of_two_equal_parts_species_in_core/CH3CF3_C2H5F/chemkin/copies/copy_chem.cti'
gas = ct.Solution(directory)
blend_name = 'CH3CF3_C2H5F'

vol_frac_list = list(np.linspace(0.025, 0.25, 30))

results = dict()


#this part calculates the flame calcs for all of the vol_frac list, using the previous vol frac count as a guess for the current one
for i in  range(len(vol_frac_list)):
    try:

        string = f'****************************starting new volume fraction **************************'
        print(string)
        
        x = vol_frac_list[i]
        norm_ox = (1-x)*.21

        vol_frac_dict = {'CH3CF3(1)': (x/2/norm_ox), 'C2H5F(188)': (x/2/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}
        gas.TPX = To, Po, vol_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        try:
            if i!=0:
                d = f'data/{vol_frac_list[i-1]}.csv'
                if os.path.exists(d):  
                    arr2 = ct.SolutionArray(gas)
                    arr2.read_csv(d)
                    flame.set_initial_guess(data=arr2)
                    print(' initial guess has been set')
        except: 
            print('initial guess not set for this volume fraction')
        
        
        #"False" stops the calculation from retrying over and over
        flame.solve(loglevel=loglevel, auto=False)
        Su = flame.velocity[0]
        results[x] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        df1.to_csv(f'data/{x}.csv', index=False)            
            
    except Exception as e: 
        print(f'********************passed volume fraction:{vol_frac_list[i]}, error: {e}*************************************')
        pass

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)

    
with open('blend_flamespeeds.csv', 'w+') as g:
    writers = csv.writer(g)
    writers.writerow(vol_fracs)
    writers.writerow(flame_speeds)
        
