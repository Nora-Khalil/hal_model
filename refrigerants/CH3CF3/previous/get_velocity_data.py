import os

directory = './data'

file_names = os.listdir(directory)

#file_names = ['.test_0.07250000000000001.csv.swp']
#file_names = ['test_0.07.csv']

#take list of file names and sort them
file_name_dict_unsorted = {}
for name in file_names:
    if name.startswith('.') and name.endswith('.swp'): 
        new_name = name.replace('.','',1).replace('.swp','') 
        file_name_dict_unsorted[new_name] = new_name[5:11]
    else: 
        file_name_dict_unsorted[name] = name[5:11]
volfrac_list_sorted = sorted(file_name_dict_unsorted.values())
  
#makes list of volume fractions
volume_fractions = []
for volfrac in volfrac_list_sorted:
    volume_fractions.append(volfrac)
    for character in volfrac:
        if character.isalpha():
            #print('there is a character')
            new_volfrac = volfrac.replace(volfrac[-2:],'')
            volume_fractions.append(new_volfrac)
            volume_fractions.remove(volfrac)


#take list of volume fractions (list of strings) and change to a list of floats
volume_fractions_float = []
for string in volume_fractions: 
    volume_fractions_float.append(float(string))
 
    
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
