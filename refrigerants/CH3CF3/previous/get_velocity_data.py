import os
import re

directory = './data'

file_names = os.listdir(directory)

file_name_dict_unsorted = {}

for name in file_names:
    match = re.match('test_([0-9.]+)\.csv', name)
    if match:
        volume_fraction = float(match.group(1))
        file_name_dict_unsorted[name] = volume_fraction
    else:
        print(f"didn't match {name}")

volume_fractions_float = sorted(file_name_dict_unsorted.values())

    
print('Volume fractions are:')
print(volume_fractions_float)


#makes list of velocities 

flamespeeds = []

filenames_sorted_dict = {k:v for k,v in sorted(file_name_dict_unsorted.items(), key=lambda item: item[1])}
for file in filenames_sorted_dict.keys():
    d = os.path.join(directory,file)
    with open(d,'r') as csvfile:
        #for characters in csvfile.readlines()[1]:
        match = re.match('0,0\.0,([0-9.]+),298\.0',csvfile.readlines()[1])
        flamespeeds.append(float(match.group(1)))
print('Flame speeds are') 
print(flamespeeds)
        