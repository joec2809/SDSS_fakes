import numpy as np
import matplotlib.pyplot as plt

arr = np.zeros((2, 10))

x_err = 3

for i in range(10):
    arr[0][i] = i
    arr[1][i] = 2*i

fig, ax = plt.subplots()

ax.errorbar(arr[0], arr[1], xerr = x_err)

plt.show()