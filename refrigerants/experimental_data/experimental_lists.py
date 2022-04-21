import os
import re
import csv

directory = '.'
d = os.path.join(directory, 'experimental_data.csv')

data = dict()

with open(d,'r') as csvfile: 
    csvreader = csv.DictReader(csvfile)

    for row in csvreader:
        formula = row['formula']
        vol_frac = row['volume_frac']
        flamespeed = row['Su (cm/s)']

        if formula not in data:
            data[formula] = ([], [])  # tuple of two lists, one for volume fractions, one for flame speeds
        data[formula][0].append(vol_frac)
        data[formula][1].append(flamespeed)
        
        print(formula, vol_frac, flamespeed)
    
print(data)


    #R-32,FCF,CH2F2,0.1736,0.9,0.15624,5.97,air,298,1,CVM,10,"T