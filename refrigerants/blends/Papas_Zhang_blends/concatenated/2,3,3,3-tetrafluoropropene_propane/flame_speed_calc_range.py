#########calculates flamespeeds for the .cti files at different equivalence ratios. Uses initial guess from the previous model #############################


#####will be used to calculate the burning velocity of ******stoichiometric******* 2,3,3,3-tetrafluoropropene-propane mix increasing mole fraction of tetra 

import cantera as ct
import numpy as np
import pandas as pd
import os
import csv 
import sys
import re


print("Running Cantera Version: " + str(ct.__version__))


To = 298
Po = ct.one_atm
directory = '/work/westgroup/nora/Code/projects/halogens/refrigerants/blends/Papas_Zhang_blends/concatenated/2,3,3,3-tetrafluoropropene_propane/copies/copy_2,3,3,3-tetrafluoropropene_propane.cti'


gas = ct.Solution(directory)

#### has tuples containing (tetra mole fraction, stoichiometric dictionary at that mole fraction)
master_dictionary = [(0.0,
  {'C3H8(165)': 0.2,
   'C3H2F4(1)': 0.0,
   'O2(3)': 1.0,
   'N2': 3.7599999999999993}),
 (0.041666666666666664,
  {'C3H8(165)': 0.19574468085106386,
   'C3H2F4(1)': 0.00851063829787234,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.08333333333333333,
  {'C3H8(165)': 0.19130434782608693,
   'C3H2F4(1)': 0.017391304347826084,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.125,
  {'C3H8(165)': 0.18666666666666668,
   'C3H2F4(1)': 0.02666666666666667,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.16666666666666666,
  {'C3H8(165)': 0.18181818181818185,
   'C3H2F4(1)': 0.03636363636363636,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.20833333333333331,
  {'C3H8(165)': 0.17674418604651163,
   'C3H2F4(1)': 0.04651162790697674,
   'O2(3)': 1.0,
   'N2': 3.7600000000000002}),
 (0.25,
  {'C3H8(165)': 0.17142857142857143,
   'C3H2F4(1)': 0.05714285714285714,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.29166666666666663,
  {'C3H8(165)': 0.16585365853658535,
   'C3H2F4(1)': 0.06829268292682925,
   'O2(3)': 1.0,
   'N2': 3.7599999999999993}),
 (0.3333333333333333,
  {'C3H8(165)': 0.16,
   'C3H2F4(1)': 0.07999999999999999,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.375,
  {'C3H8(165)': 0.15384615384615385,
   'C3H2F4(1)': 0.09230769230769231,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.41666666666666663,
  {'C3H8(165)': 0.1473684210526316,
   'C3H2F4(1)': 0.10526315789473682,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.4583333333333333,
  {'C3H8(165)': 0.14054054054054055,
   'C3H2F4(1)': 0.11891891891891891,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.5,
  {'C3H8(165)': 0.13333333333333333,
   'C3H2F4(1)': 0.13333333333333333,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.5416666666666666,
  {'C3H8(165)': 0.12571428571428572,
   'C3H2F4(1)': 0.14857142857142855,
   'O2(3)': 1.0,
   'N2': 3.7600000000000002}),
 (0.5833333333333333,
  {'C3H8(165)': 0.11764705882352942,
   'C3H2F4(1)': 0.16470588235294115,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.625,
  {'C3H8(165)': 0.10909090909090909,
   'C3H2F4(1)': 0.18181818181818182,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.6666666666666666,
  {'C3H8(165)': 0.1,
   'C3H2F4(1)': 0.19999999999999998,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.7083333333333333,
  {'C3H8(165)': 0.0903225806451613,
   'C3H2F4(1)': 0.21935483870967737,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.75, {'C3H8(165)': 0.08, 'C3H2F4(1)': 0.24, 'O2(3)': 1.0, 'N2': 3.76}),
 (0.7916666666666666,
  {'C3H8(165)': 0.06896551724137932,
   'C3H2F4(1)': 0.26206896551724135,
   'O2(3)': 1.0,
   'N2': 3.7599999999999993}),
 (0.8333333333333333,
  {'C3H8(165)': 0.05714285714285716,
   'C3H2F4(1)': 0.28571428571428564,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.875,
  {'C3H8(165)': 0.044444444444444446,
   'C3H2F4(1)': 0.3111111111111111,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.9166666666666666,
  {'C3H8(165)': 0.03076923076923078,
   'C3H2F4(1)': 0.3384615384615384,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.9583333333333333,
  {'C3H8(165)': 0.016000000000000028,
   'C3H2F4(1)': 0.36799999999999994,
   'O2(3)': 1.0,
   'N2': 3.7600000000000002}),
 (1.0,
  {'C3H8(165)': 0.0,
   'C3H2F4(1)': 0.4,
   'O2(3)': 1.0,
   'N2': 3.7599999999999993})]

results = {}

for fraction, dictionary in master_dictionary:
    try: 
        
        tol_ss = [1.0e-13, 1.0e-9]  #abs and rel tolerances for steady state problem
        tol_ts = [1.0e-13, 1.0e-9]  #abs and rel tie tolernces for time step function
        
       
        print(f'****************************starting new mole fraction for 2,3,3,3-tetrafluoropropene: {fraction}**************************')

        mole_frac_dict = dictionary

        print(mole_frac_dict)
        gas.TPX = To, Po, mole_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.flame.set_steady_tolerances(default=tol_ss)   #set tolerances
        flame.flame.set_transient_tolerances(default=tol_ts)
        #flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.set_refine_criteria(ratio=5, slope=0.25, curve=0.27)
        flame.max_time_step_count = 900
        loglevel = 1 
        
        ######################################################################################
#         if i!=0:
#             d = f'./data/range10pts_test_{BTPmole_list[i-1]}.csv'
#             if os.path.exists(d):  
#                 arr2 = ct.SolutionArray(gas)
#                 arr2.read_csv(d)
#                 flame.set_initial_guess(data=arr2)
#                 print(' initial guess has been set')
        #######################################################################################        
                
        #"False" stops the calculation from retrying over and over, thanks Chao 
        #flame.solve(loglevel=loglevel, auto=False)
        flame.solve(loglevel=loglevel, auto=True)
        Su = flame.velocity[0]

        results[fraction] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        #edited this here!! index=False
        df1.to_csv(f'{fraction}.csv', index=False)
    except Exception as e: 
        print(f'********************   passed mole fraction:{fraction}, error: {e}    *************************************')
        pass


mole_fraction_of_tetra = list(results.keys())
flame_speeds = list(results.values())


print("mole fractions of 2,3,3,3-tetrafluoropropene are:")
print(mole_fraction_of_tetra)

print("flame speeds are:")
print(flame_speeds)


with open(f'flame_speeds.csv', 'w+') as g:
    g.write(directory)
    g.write('\n')
    writers = csv.writer(g)
    writers.writerow(mole_fraction_of_tetra)
    writers.writerow(flame_speeds)

        

