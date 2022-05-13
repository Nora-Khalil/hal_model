import os 
import re

directory = '../'

list_of_blends = [file for file in os.listdir(directory) if re.match('[A-Z0-9]+_[A-Z0-9]+', file) ]

i=1

#delete duplicates in list_of_blends
for idx,x in enumerate(list_of_blends): 
    match1 = re.match('([A-Z0-9]+)_([A-Z0-9]+)', x)
    other_blends = list_of_blends[:idx] + list_of_blends[idx+1:]
    for idc, c in enumerate(other_blends): 
        match2 = re.match('([A-Z0-9]+)_([A-Z0-9]+)', c)
        if (match1.group(1)==match2.group(2)) and (match1.group(2)==match2.group(1)):
            print(match1.group(1), match2.group(2), idc)
            list_of_blends.remove(c)
            i +=1
            print(i)


#create list of files that will need to be deleted
list_of_blends_initially = [file for file in os.listdir(directory) if re.match('[A-Z0-9]+_[A-Z0-9]+', file) ]
to_be_deleted = [x for x in list_of_blends_initially if x not in list_of_blends] #remember, list of blends was updated at this point to remove all duplicates

#delete using terminal commands
for i in to_be_deleted: 
    os_command1 = f'rm -rf ../{i}'
    os.system(os_command1)

    
list_of_blends_after = [file for file in os.listdir(directory) if re.match('[A-Z0-9]+_[A-Z0-9]+', file) ]   

#should get list of 21 files
print(len(list_of_blends_after))
