import os

directory = '../'

species = {'C2H6': [0.025,0.1525,0.0566],'CH4': [0.05, 0.1525,0.095] ,'CH2F2': [0.11,0.225,0.174] ,'CH3F': [0.08,.225,0.123],'C2H5F': [0.025,0.150,0.0654],'CH2FCH2F' :[0.0375,0.150,0.0775],'CH2FCHF2' :[0.075,0.150,0.095],'CH3CF3' :[0.07,0.150,0.095],'CH3CHF2':[0.0375, 0.150,0.0775]}


flame_speed_calc_file = '''

import cantera as ct
import numpy as np
import pandas as pd

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

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
    mole_frac_dict = {'***label***': (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
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
    pd.to_csv(f'data/test_{x}.csv')

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)

print(xresults)
print("volume fractions list is " + str(mole_frac_list))

'''.replace('***label***',label) 
   return block

def make_bash_script(label):
    bash_script = """#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=unlimited
#SBATCH --job-name=fc_{label}
#SBATCH --error=fc.slurm.log
#SBATCH --output=fc_output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc.py
""".replace('label',label)

    return bash_script

if __name__ == "__main__":

    for label,values in species.items():  
       d = os.path.join(directory,label,'flame_speed_calc.py')
       if os.path.exists(d):
          continue
       flame_speed_calc = flame_speed_calc_file
       flame_speed_calc += make_vol_frac_list(values[0],values[1],values[2])
       flame_speed_calc +=make_for_loop(label)
       bash_script = make_bash_script(label)
       d = os.path.join(directory,label)
       os.makedirs(d,exist_ok=True)
       with open(os.path.join(d,'flame_speed_calc.py'),'w') as f: 
         for l in flame_speed_calc: 
            f.write(l)
       with open(os.path.join(d,'flame_speed_run.sh'),'w') as f: 
         for l in bash_script: 
            f.write(l)
       par_dir= os.path.join(directory,label)
       if os.path.exists(par_dir): 
          continue 
       end_dir = 'data'
       data = os.path.join(par_dir,end_dir,'')
       os.makedirs(data)  
