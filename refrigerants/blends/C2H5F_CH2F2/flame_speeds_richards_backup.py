import cantera as ct
import numpy as np
import pandas as pd
import os 


print("Running Cantera Version: " + str(ct.__version__))



To = 298
Po = 1e5 # ct.one_atm

gas = ct.Solution('./chemkin/chem_new.cti')


mole_frac_list = np.arange(0.025, 0.25, step=0.025)
#mole_frac_list = [0.1, 0.125, 0.15, 0.175]
#mole_frac_list = [0.15]



#mole_frac_list = np.arange(0.025, 0.25, step=0.025)
print(mole_frac_list)


results = {}

#this part does the rest of the list, using the previous vol frac count as a quess for the current one
for i in  range(len(mole_frac_list)):
    try:
        # x is mole fraction of fuel, which is 50% each component
        print(f'************************starting calculation for volume fraction:{mole_frac_list[i]}*****************')
        x = mole_frac_list[i]
        norm_ox = (1-x)*.21
        mole_frac_dict = {'C2H5F(1)': (x/2/norm_ox), 
                          'CH2F2(2)': (x/2/norm_ox), 
                          'O2(3)':((1-x)*.21)/norm_ox, 
                          'N2':((1-x)*0.79)/norm_ox } 
        gas.TPX = To, Po, mole_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        if i!=0:
          d = f'data_richards/test_{mole_frac_list[i-1]}.csv'
          if os.path.exists(d):  
            arr2 = ct.SolutionArray(gas)
            arr2.read_csv(d)
            flame.set_initial_guess(data=arr2)
        flame.solve(loglevel=loglevel, auto=False)
        #flame.solve(loglevel=loglevel, auto=True)
        Su = flame.velocity[0]
        results[x] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        df1.to_csv(f'data_richards/test_{x}.csv')

         #delete first column of csv file to avoid error
        df2 = pd.read_csv(f'data_richards/test_{x}.csv')
        # If you know the name of the column skip this
        first_column = df2.columns[0]
        # Delete first
        df2 = df2.drop([first_column], axis=1)
        df2.to_csv(f'data_richards/test_{x}.csv', index=False)
    except Exception as e: 
        print(f'********************passed volume fraction:{mole_frac_list[i]}, error: {e}*************************************')
        pass

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)



