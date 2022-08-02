
        ############################# calculates flamespeeds for the .cti files at different equivalence ratios. Uses initial guess from the previous model #############################

import cantera as ct
import numpy as np
import pandas as pd
import os
import csv 

print("Running Cantera Version: " + str(ct.__version__))


To = 298
Po = 1e5 # ct.one_atm
vol_frac_list = np.arange(0.025, 0.25, step=0.01)

cti_files = ['dup_chem.cti']

for file in cti_files: 
    #make directory to store flame speed calculations, and a header for csv files
    d = f'./flame_calcs/chem_CALC'
    os.makedirs(d, exist_ok=True)
    gas = ct.Solution(f'./chemkin/dups/{file}')
    results = {}    
    #this part calculates the flame calcs for all of the vol_frac list, using the previous vol frac count as a guess for the current one
    for i in  range(len(vol_frac_list)):
        try:

            string = f'****************************starting new volume fraction: {vol_frac_list[i]}**************************'
            print(string)

            x = vol_frac_list[i]
            norm_ox = (1-x)*.21

        
            vol_frac_dict = {'C2H4F2(1)': (x/2/norm_ox), 'CH3F(2)': (x/2/norm_ox), 'O2(3)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}
            gas.TPX = To, Po, vol_frac_dict
            width = 0.08
            flame = ct.FreeFlame(gas, width=width)
            flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
            flame.max_time_step_count = 900
            loglevel = 1 
            if i!=0:
                d = f'./flame_calcs/chem_CALC/test_{vol_frac_list[i-1]}.csv'
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
            #edited this here!! index=False
            df1.to_csv(f'./flame_calcs/chem_CALC/test_{x}.csv', index=False)


        except Exception as e: 
            print(f'********************passed volume fraction:{vol_frac_list[i]}, error: {e}*************************************')
            pass

    vol_fracs = list(results.keys())
    flame_speeds = list(results.values())


    print("volume fractions are:")
    print(vol_fracs)

    print("flame speeds are:")
    print(flame_speeds)


    with open('final_calcs.csv', 'w+') as g:
        writers = csv.writer(g)
        writers.writerow(vol_fracs)
        writers.writerow(flame_speeds)


        