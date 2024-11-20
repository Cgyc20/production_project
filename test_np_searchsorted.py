import numpy as np


xarray = np.arange(0,10)*0.1

x = 0.01
print(xarray)

print(np.searchsorted(xarray,x))