

import sys
sys.path.insert(0, '/home/khalil.nor/cantera/build/python')
import cantera as ct
import os


print("Running Cantera Version: " + str(ct.__version__))

#Constants
To = 298
Po = ct.one_atm

gas = ct.Solution('./cantera/chem.cti')


i = 0.075
mole_frac_list = []

while i < 0.15: 
   mole_frac_list.append(i)
   i += 0.0025 

#can use below to only test one vol_frac
#mole_frac_list = [0.095]


results = {}

for x in mole_frac_list: 
    norm_ox = (1-x)*.21
    mole_frac_dict = {'C2H3F3(1)': (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
    gas.TPX = To, Po, mole_frac_dict
    width = 0.08
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
    flame.max_time_step_count = 900
    loglevel = 1 
    flame.transport_model='Multi'
    flame.solve(loglevel=loglevel, auto=True)
    Su = flame.velocity[0]
    results[x] = Su
    sltn = flame.to_solution_array()
    
    #save to csv file
    directory = './data_tp_multi'
    with open(os.path.join(directory,f'{x}.csv'),'w') as f:
       f.write(str(Su))
       f.write(str(results.keys()))
       f.write(str(results.values()))
    

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)

#Test for different transport properties
properties = ['binary_diff_coeffs', 'electrical_conductivity', 'get_binary_diff_coeffs_polynomial(1,2)', 'get_collision_integral_polynomials(1,2)', 'get_thermal_conductivity_polynomial(1)', 'get_viscosity_polynomial(1)', 'mix_diff_coeffs', 'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole', 'mobilities', 'multi_diff_coeffs', 'species_viscosities', 'thermal_conductivity', 'thermal_diff_coeffs']

for prop in properties: 
   string = f'gas.{prop}'
   try: 
       print(prop)
       eval(string)
   except AttributeError: 
      pass 
      print(f'gas object has no {prop} attribute')
   except NotImplementedError: 
      pass
      print(f'gas object has no {prop} attribute')
   finally: 
      continue 


