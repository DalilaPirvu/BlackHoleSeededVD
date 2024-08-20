# To run script type:
#### python3 script.py >> output.txt ####
# into terminal

import numpy as np
from numpy.linalg import eigh
from findSpectrum import *

frac = 2.

for lamb in [0.028]:
    for nLat in 2**np.arange(9,10):
        for lenLat in [100]:
 
            spectrum = Spectrum(nLat, lenLat, lamb, frac)
            spectrum.writeFrequencies()

print('End of session.')
print('/n')
        
