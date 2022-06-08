import os

directory = '../'

species = {'C2H6': [0.025,0.1525,0.0566],'CH4': [0.05, 0.1525,0.095] ,'CH2F2': [0.11,0.225,0.174] ,'CH3F': [0.08,.225,0.123],'C2H5F': [0.025,0.150,0.0654],'CH2FCH2F' :[0.0375,0.150,0.0775],'CH2FCHF2' :[0.075,0.150,0.095],'CH3CF3' :[0.07,0.150,0.095],'CH3CHF2':[0.0375, 0.150,0.0775]}
#species = {'C2H6': [0.025,0.1525,0.0566]}

flame_speed_calc_file = '''

import sys
sys.path.insert(0, '/home/khalil.nor/cantera/build/python')
import cantera as ct
import os


print("Running Cantera Version: " + str(ct.__version__))

#Constants
To = 298
Po = ct.one_atm

gas = ct.Solution('./cantera/chem.cti')

'''

def make_vol_frac_list(vol_frac_lower_bound, vol_frac_upper_bound,equiv_rat_1):

    block = \
   '''
i = ***vol_frac_lower_bound***
mole_frac_list = []

while i < ***vol_frac_upper_bound***: 
   mole_frac_list.append(i)
   i += 0.0025 

#can use below to only test one vol_frac
#mole_frac_list = [***equiv_rat_1***]

'''.replace('***vol_frac_lower_bound***',str(vol_frac_lower_bound)).replace('***vol_frac_upper_bound***',str(vol_frac_upper_bound)).replace('***equiv_rat_1***',str(equiv_rat_1))
    
    return block

def make_for_loop(label):

   block = \
   '''
results = {}

for x in mole_frac_list: 
    norm_ox = (1-x)*.21
    mole_frac_dict = {'***label***(1)': (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
    gas.TPX = To, Po, mole_frac_dict
    width = 0.08
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
    flame.max_time_step_count = 900
    loglevel = 1 
    flame.transport_model='Mix'
    flame.solve(loglevel=loglevel, auto=True)
    Su = flame.velocity[0]
    results[x] = Su
    sltn = flame.to_solution_array()
    
    #save to csv file
    directory = './data_tp'
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


'''.replace('***label***',label) 
   return block

def make_bash_script(label):
    bash_script = """#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=12-00:00:00
#SBATCH --job-name={label}_mix
#SBATCH --error=tp_mix.slurm.log
#SBATCH --output=tp_mix.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate proj_transport
PYTHONPATH=/home/khalil.nor/cantera/build/python python tp_mix_flame_speed_calc.py
##python tp_mix_flame_speed_calc.py
""".replace('label',label)

    return bash_script

if __name__ == "__main__":

    for label,values in species.items():  
       d = os.path.join(directory,label,'tp_mix_flame_speed_calc.py')
       if os.path.exists(d):
          continue
       flame_speed_calc = flame_speed_calc_file
       flame_speed_calc += make_vol_frac_list(values[0],values[1],values[2])
       flame_speed_calc +=make_for_loop(label)
       bash_script = make_bash_script(label)
       d = os.path.join(directory,label)
       with open(os.path.join(d,'tp_mix_flame_speed_calc.py'),'w') as f: 
         for l in flame_speed_calc: 
            f.write(l)
       with open(os.path.join(d,'run_tp_mix_flame_speed_calc.sh'),'w') as f: 
         for l in bash_script: 
            f.write(l)
       data_tp = f'../{label}/data_tp'
       os.makedirs(data_tp,exist_ok=True)
       """ par_dir= os.path.join(directory,label)
       if os.path.exists(par_dir): 
          continue 
       end_dir = 'data_tp'
       data = os.path.join(par_dir,end_dir,'')
       os.makedirs(data)   """
