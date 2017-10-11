import sys
import numpy as np

filename = sys.argv[1]

vals = np.loadtxt(filename,skiprows=1,dtype=float,delimiter=',')

for i,v in enumerate(vals):

    if i>0:
        if v[0]==vals[i-1][0]:
            vals[i-1][0] += 0.5
            v[0]         -= 0.5

            #print vals

for v in vals:
    output = ','.join(map(str,v)) 
    print(output)




    



