import sys
sys.path.insert(0, '/home/khalil.nor/cantera/build/python')
import cantera as ct 
import os

#create gas object
gas = ct.Solution('./cantera/chem.cti')

To = 298
Po = ct.one_atm

#create flame object, can evaluate transport properties at just one equivalence ratio or many 

volume_fraction_list = [0.095] #test only when volume fraction is 9.5% and equivalence ratio = 1

#use while loop below to determine properties across all volume fractions
'''
i = 0.05 #lower bound of volume_fraction
volume_fraction_list = []

while i < 0.1525:
   volume_fraction_list.append(i)
   i += 0.0025

'''

results = {}

#Test for different transport properties
properties = ['binary_diff_coeffs', 'electrical_conductivity', 'get_binary_diff_coeffs_polynomial(1,2)', 'get_collision_integral_polynomials(1,2)', 'get_thermal_conductivity_polynomial(1)', 'get_viscosity_polynomial(1)', 'mix_diff_coeffs', 'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole', 'mobilities', 'multi_diff_coeffs', 'species_viscosities', 'thermal_conductivity', 'thermal_diff_coeffs']


for vol in volume_fraction_list:
    norm_ox = (1-vol)*.21
    equivalence_ratio_dict = {'CH4(1)': (vol/norm_ox), 'O2(2)':((1-vol)*.21)/norm_ox, 'N2':((1-vol)*0.79)/norm_ox }
    gas.TPX = To, Po, equivalence_ratio_dict
    width = 0.08
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.max_time_step_count = 900
    loglevel = 1
    
#choose transport model to use, either "Mix" or "Multi"
    flame.transport_model='Mix'
    flame.solve(loglevel=loglevel, auto=True)

#now test for properties in gas flame
    print('TRANSPORT PROPERTIES OF GAS')
    for prop in properties: 
       try: 
           print(prop)
           eval(f'print(gas.{prop})')
           print()
       except AttributeError: 
          print(f'gas object has no {prop} attribute')
          print()
       except NotImplementedError: 
          print(f'gas object has no {prop} attribute')
          print()  
       except ct.CanteraError as err:
          print(f'Cantera error for {prop} attribute: {err}')
          print()
        
#now test for properties in free flame
    print('TRANSPORT PROPERTIES OF FREE FLAME')
    for prop in properties: 
       try: 
           print(prop)
           eval(f'print(gas.{prop})')
           print()
       except AttributeError: 
          print(f'gas object has no {prop} attribute')
          print()
       except NotImplementedError: 
          print(f'gas object has no {prop} attribute')
          print()  
       except ct.CanteraError as err:
          print(f'Cantera error for {prop} attribute: {err}')
          print()
        

