import os

directory = '../'

species = ['C2H5F','C2H6','CH2F2','CH2FCH2F','CH2FCHF2','CH3CF3','CH3CHF2','CH3F','CH4']

def make_bash_script(label): 
 bash_script = """#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --job-name=gd_label
#SBATCH --error=gd.slurm.log
#SBATCH --output=gd_output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west



python get_velocities_from_data.py
""".replace('label',label)

    return bash_script

get_velocity_from_data = '''

import os




with open(./test




'''




if __name__ == "__main__":

    for label in species:
       #make get_velocities_fie in labels folder 
       label_directory = os.path.join(directory,label)
       get_velocities_file = os.path.join(label_directory,'data','get_velocities.py')
       os.makedirs(get_velocities_file)

       #make get_velocities_fie in previous folder 
       previous_folder = os.path.join(label_directly,'previous')
       get_velocities_file_in_previous = os.path.join(previous_folder,'data','get_velocities.py')
       os.makedirs(get_velocities_file_in_previous)
       
       #write to get_velocities_file in label
       with open(get_velocities_file,'w') as f: 
         for l in get_velocity_from_data: 
           f.write(l)

       #write to get_velocities_file in previous
       with open(get_velocities_file_in_previous,'w') as f: 
         for l in get_velocity_from_data: 
           f.write(l)

       #make bash scripts in label and in previous
       bash_script = make_bash_script(label)
       with open(os.path.join(label_directory,'data','get_velocities_run.sh'),'w') as f: 
         for l in bash_script: 
           f.write(l)
      with open(os.path.join(previous_folder,'data','get_velocities_run.sh'),'w') as f:
         for l in bash_script:
           f.write(l)
       
       




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
       data = os.path.join(par_di
