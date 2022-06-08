from os import uname
import cantera as ct
import numpy as np
import pandas as pd
import os

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

gas = ct.Solution('./chemkin/dups/dup_chem0076.cti')



vol_frac_list = np.arange(0.025, 0.25, step=0.025)

results = {}

for i in  range(len(vol_frac_list)):
    try:

        string = f'****************************starting new volume fraction: {vol_frac_list[i]}**************************'
        print(string)

        x = vol_frac_list[i]
        norm_ox = (1-x)*.21
        vol_frac_dict = {'C2H5F(1)': (x/2/norm_ox), 
                          'CH2F2(2)': (x/2/norm_ox), 
                          'O2(3)':((1-x)*.21)/norm_ox, 
                          'N2':((1-x)*0.79)/norm_ox } 
        gas.TPX = To, Po, vol_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        try:
            if i!=0:
                d = f'data_all/test_{vol_frac_list[i-1]}.csv'
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
        #flame.solve(loglevel=loglevel, auto=True)
        Su = flame.velocity[0]
        results[x] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        df1.to_csv(f'data_all/test_{x}.csv')



         #delete first column of csv file to avoid error
        df2 = pd.read_csv(f'data_all/test_{x}.csv')
        first_column = df2.columns[0]
        # Delete first
        df2 = df2.drop([first_column], axis=1)
        df2.to_csv(f'data_all/test_{x}.csv', index=False)
    except Exception as e: 
        print(f'********************passed volume fraction:{vol_frac_list[i]}, error: {e}*************************************')
        pass

    vol_fracs = list(results.keys())
    flame_speeds = list(results.values())


    print("volume fractions are:")
    print(vol_fracs)
 
    print("flame speeds are:")
    print(flame_speeds)
    

