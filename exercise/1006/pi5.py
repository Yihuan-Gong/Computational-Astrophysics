import numpy as np

# Nunber of sample points
N = 10000

# Random seed for x and y
seedx = 100
seedy = 200

# Gnerate independent random series for x and y
np.random.seed(seedx)
x = np.random.random(N)

np.random.seed(seedy)
y = np.random.random(N)

# Count points within the circle
r = np.sqrt(x**2 + y**2)
in_circle = len(r[r<1])

# Area of the circle
A = 4*in_circle/N

print('Area of unit circle = ' + str(A))

