import os
import re

directory = '.'
d= os.path.join(directory,'experimental_data.csv')

with open(d,'r') as csvfile: 
    lines = csvfile.readlines()
    for line in lines[1:]:
        match = re.match('R-[a-zA-Z0-9]+,[^,]+,([\w]+)[\s+]?,+[0-9.]+,[0-9.]+,([0-9.]+),([0-9.]+),',line)
        #match = re.match('R-[a-zA-Z0-9]+,[A-Z()=]+,([\w]+),[0-9.]+,[0-9.]+,([0-9.]+),([0-9.]+),',line)
        if match: 
            group, vol_frac, flamespeed = (match.group(1)), float(match.group(2)), float(match.group(3))
            print(group, vol_frac, flamespeed)
        else: 
            print('code better Nora')
 






    #R-32,FCF,CH2F2,0.1736,0.9,0.15624,5.97,air,298,1,CVM,10,"T