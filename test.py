import numpy as np

from astropy import units as u

arr = np.arange(1,10,1)*u.Unit('m s-1')

print(type(arr.value))