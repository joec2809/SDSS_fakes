import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt

with fits.open("./dr16/sdss/spectro/redux/26/spectra/0854/spec-0854-52373-0369-0.9.fits", comments='#') as hdul:
            header = hdul[0].header
            spectrum_data = hdul[1].data

#Assumes the spectrum data is in the format:
#Wavelength, Flux, FluxError

wave = np.array(spectrum_data['loglam'])
flux = np.array(spectrum_data['flux'])

fig, ax = plt.subplots()
ax.plot(wave,flux)

plt.show()