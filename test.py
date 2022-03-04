import numpy as np

arr = np.zeros((2, 10))

for i in range(10):
    arr[0][i] = i
    arr[1][i] = 2*i

print(arr[:,0])