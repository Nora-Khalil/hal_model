

import cantera as ct
import numpy as np
import pandas as pd

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

gas = ct.Solution('./chemkin/chem_final.cti')


i = 0.08
mole_frac_list = []

while i < 0.15: 
   mole_frac_list.append(i)
   i += 0.0025 



results = {}

for x in mole_frac_list: 
    norm_ox = (1-x)*.21
    mole_frac_dict = {'C2H3F3(1)': (x/norm_ox), 'CH3F(2)': (x/norm_ox), 'O2(3)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
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
    pd.to_csv(f'data_6d/test_{x}.csv')

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)

