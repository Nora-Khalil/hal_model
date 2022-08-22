#########calculates flamespeeds for the .cti files at different equivalence ratios. Uses initial guess from the previous model #############################

import cantera as ct
import numpy as np
import pandas as pd
import os
import csv 
print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = ct.one_atm
directory = '/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/cantera/Nora/2_BTP/FFCM_seed/CH4_with_seed/chemkin/copies/copy_chem.cti'
gas = ct.Solution(directory)


#vol_frac_list = np.arange(0.5, 1.2, step=0.07)


#can use below to only test one vol_frac
vol_frac_list = [0.095]


results = {}

for i in  range(len(vol_frac_list)):
    try: 
        
        tol_ss = [1.0e-13, 1.0e-9]  #abs and rel tolerances for steady state problem
        tol_ts = [1.0e-13, 1.0e-9]  #abs and rel tie tolernces for time step function
        
        x = vol_frac_list[i]
        norm_ox = (1-x)*.21
        
        
        print(f'****************************starting new volume fraction: {x}**************************')

        vol_frac_dict = {'CH4(1)': (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}
        print(vol_frac_dict)
        print(f"O2/CH4 ratio = {vol_frac_dict['O2(2)']/vol_frac_dict['CH4(1)']}. Complete combustion takes 2")
        gas.TPX = To, Po, vol_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.flame.set_steady_tolerances(default=tol_ss)   #set tolerances
        flame.flame.set_transient_tolerances(default=tol_ts)
        #flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.set_refine_criteria(ratio=5, slope=0.25, curve=0.27)
        flame.max_time_step_count = 900
        loglevel = 1 
#         if i!=0:
#             d = f'./data/David_test_{vol_frac_list[i-1]}.csv'
#             if os.path.exists(d):  
#                 arr2 = ct.SolutionArray(gas)
#                 arr2.read_csv(d)
#                 flame.set_initial_guess(data=arr2)
#                 print(' initial guess has been set')
        #"False" stops the calculation from retrying over and over, thanks Chao 
        #flame.solve(loglevel=loglevel, auto=False)
        flame.solve(loglevel=loglevel, auto=True)
        Su = flame.velocity[0]
        results[x] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        #edited this here!! index=False
        df1.to_csv(f'./data/CH4_with_FFCM_seed_{x}.csv', index=False)
    except Exception as e: 
        print(f'********************passed volume fraction:{vol_frac_list[i]}, error: {e}*************************************')
        pass


vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)


with open('final_calcs_CH4_with_FFCM_seed.csv', 'w+') as g:
    g.write(directory)
    g.write('\n')
    writers = csv.writer(g)
    writers.writerow(vol_fracs)
    writers.writerow(flame_speeds)

        

