
    
#########calculates flamespeeds for the .cti files at different equivalence ratios. Uses initial guess from the previous model #############################

from os import uname
import cantera as ct
import numpy as np
import pandas as pd
import os
import csv 

print("Running Cantera Version: " + str(ct.__version__))

    
To = 298
Po = 1e5 # ct.one_atm

vol_frac_list = list(np.linspace(0.025, 0.25, 30))

gas = ct.Solution(f'./cantera/chem.cti')

results = {}  

  
#this part calculates the flame calcs for all of the vol_frac list, using the previous vol frac count as a guess for the current one
for i in  range(len(vol_frac_list)):
    try:
        
        x = vol_frac_list[i]
        norm_ox = (1-x)*.21
        
        print(f'****************************starting new volume fraction: {x}**************************')
        

        vol_frac_dict = {'C2H5F(1)': (x/2/norm_ox), 'CH2F2(2)': (x/2/norm_ox), 'O2(3)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}
    
        print(vol_frac_dict)
        
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
        except Exception as e:
            print(f'******************************will try to solve by scratch because of error: {e}******************************************')
            pass
        #"False" stops the calculation from retrying over and over, thanks Chao 
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
    
with open('final_calcs.csv', 'w+') as g:
    writers = csv.writer(g)
    writers.writerow(vol_fracs)
    writers.writerow(flame_speeds)


To = 298
Po = 1e5 # ct.one_atm

vol_frac_list = list(np.linspace(0.025, 0.25, 30))

gas = ct.Solution(f'./cantera/chem.cti')

results = {}  

  
#this part calculates the flame calcs for all of the vol_frac list, using the previous vol frac count as a guess for the current one
for i in  range(len(vol_frac_list)):
    try:
        
        x = vol_frac_list[i]
        norm_ox = (1-x)*.21
        
        print(f'****************************starting new volume fraction: {x}**************************')
        

        vol_frac_dict = {'C2H5F(1)': (x/2/norm_ox), 'C2H4F2(2)': (x/2/norm_ox), 'O2(3)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}
    
        print(vol_frac_dict)
        
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
        except Exception as e:
            print(f'******************************will try to solve by scratch because of error: {e}******************************************')
            pass
        #"False" stops the calculation from retrying over and over, thanks Chao 
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
    
with open('final_calcs.csv', 'w+') as g:
    writers = csv.writer(g)
    writers.writerow(vol_fracs)
    writers.writerow(flame_speeds)


To = 298
Po = 1e5 # ct.one_atm

vol_frac_list = list(np.linspace(0.025, 0.25, 30))

gas = ct.Solution(f'./cantera/chem.cti')

results = {}  

  
#this part calculates the flame calcs for all of the vol_frac list, using the previous vol frac count as a guess for the current one
for i in  range(len(vol_frac_list)):
    try:
        
        x = vol_frac_list[i]
        norm_ox = (1-x)*.21
        
        print(f'****************************starting new volume fraction: {x}**************************')
        

        vol_frac_dict = {'C2H5F(1)': (x/2/norm_ox), 'C2H3F3(2)': (x/2/norm_ox), 'O2(3)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox}
    
        print(vol_frac_dict)
        
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
        except Exception as e:
            print(f'******************************will try to solve by scratch because of error: {e}******************************************')
            pass
        #"False" stops the calculation from retrying over and over, thanks Chao 
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
    
with open('final_calcs.csv', 'w+') as g:
    writers = csv.writer(g)
    writers.writerow(vol_fracs)
    writers.writerow(flame_speeds)

