#!/usr/bin/env python
# coding: utf-8

"""
Runs a single flame speed and saves the answer in speeds directory.
Set SLURM_ARRAY_TASK_ID to specify the species index that is modified.


"""

import cantera as ct
import numpy as np
import os

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm
volume_fraction_fuel = 0.123

gas = ct.Solution('./cantera/chem.cti')

# Modify the species
index = int(os.getenv('SLURM_ARRAY_TASK_ID', default='0'))


s = gas.species(index)
old = s.transport.well_depth
new = 2 * old
print(f"Changing species {index} {gas.species_name(index)} well depth from {old} to {new} J")
s.transport.well_depth = new


volume_fraction_oxygen = (1-volume_fraction_fuel)*.21

mole_frac_dict = {'CH3F(1)': volume_fraction_fuel,
                  'O2(2)': volume_fraction_oxygen,
                  'N2': (1 - volume_fraction_fuel - volume_fraction_oxygen) } 
gas.TPX = To, Po, mole_frac_dict
width = 0.08 # m
flame = ct.FreeFlame(gas, width=width)
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
flame.max_time_step_count = 900
loglevel = 1
flame.transport_model = 'Mix' # or 'Multi'
print(f"Using {flame.transport_model} transport model")


flame.solve(loglevel=loglevel, auto=True)

flame_speed = flame.velocity[0]

os.makedirs('speeds', exist_ok=True)

with open(f'speeds/speed{index:03}.txt','w') as f:
    f.write(f"{flame_speed}/n{gas.species_name(index)}\n")

