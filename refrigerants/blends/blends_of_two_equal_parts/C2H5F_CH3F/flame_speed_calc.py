

import cantera as ct
import numpy as np
import pandas as pd

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

#gas = ct.Solution('./cantera/chem.cti')
gas = ct.Solution('./chemkin/chemcopy.cti')

mole_frac_list = np.arange(0.025, 0.25, step=0.025)

results = {}

for x in mole_frac_list: 
    try:
        norm_ox = (1-x)*.21
        mole_frac_dict = {'C2H5F(1)': (x/2/norm_ox), 'CH3F(2)': (x/2/norm_ox), 'O2(3)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
        gas.TPX = To, Po, mole_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        flame.solve(loglevel=loglevel, auto=True)
        Su = flame.u[0]
        results[x] = Su
        sltn = flame.to_solution_array()
        pd = sltn.to_pandas()
        pd.to_csv(f'data_all/test_{x}.csv')
    except: 
        pass

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)

