import numpy as np

sum = 0
amin = 1
amax = 100
arr = np.arange(amin, amax+1, 1)
for i in arr:
    sum += i

print ('Summation from '+ str(float(amin)) + ' to ' + str(float(amax)) + ' = ' + str(float(sum)))