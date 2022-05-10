import os
import re

directory = './data'

file_names = os.listdir(directory)

#file_names = ['.test_0.07250000000000001.csv.swp']
#file_names = ['test_0.07.csv']

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
        secondword = [] #will only use the second velocity in this list 
        for line in csvfile: 
            secondword.append(line.split(',',3)[2])
        flamespeeds.append(secondword[1])


#take list of flame speeds (list of strings) and change to a list of floats
flamespeeds_float = []
for string in flamespeeds:
    flamespeeds_float.append(float(string))
  
print('Flame speeds are') 

print(flamespeeds_float)  
