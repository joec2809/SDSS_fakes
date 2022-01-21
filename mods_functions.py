import sys
import os
import math

import numpy as np
import scipy.constants as constants
import Hirogen_Functions
from astropy.io import fits
from astropy import units as u

c = constants.value('speed of light in vacuum') / 1000

def peak_selector(flux, wavelengths, line_names):
    object_name = ''
    flux_without_peaks = np.zeros(len(flux))
    for i, value in enumerate(flux):
        flux_without_peaks[i] = value
    lines = Hirogen_Functions.lines_for_analysis()
    for ii, item in enumerate(lines):
        if item in line_names:
            line_location = lines[item][0]
            shift = (((np.array(wavelengths) * c) / line_location) - c)
            shift_region, wave_region, flux_region, error_region = Hirogen_Functions.region_cutter(
                    shift, 
                    wavelengths,
                    flux*u.Unit('erg cm-2 s-1 AA-1'),
                    line_location - 40,
                    line_location + 40,
                    mode='Wavelength'
                )
            continuum = Hirogen_Functions.continua_maker(
                    wave_region,
                    flux_region,
                    # shift=Shift_Region,
                    item,
                    line_location,
                    object_name
                )
            peak_start = find_nearest(wavelengths, continuum[2][0])
            flux_without_peaks[peak_start:peak_start+len(continuum[0])] = continuum[0]

    flux -= flux_without_peaks
    return flux

def peak_remover(flux, wavelengths, line_names):
    object_name = ''
    flux_without_peaks = np.zeros(len(flux))
    for i, value in enumerate(flux):
        flux_without_peaks[i] = value
    lines = Hirogen_Functions.lines_for_analysis()
    for ii, item in enumerate(lines):
        if item in line_names:
            line_location = lines[item][0]
            shift = (((np.array(wavelengths) * c) / line_location) - c)
            shift_region, wave_region, flux_region, error_region = Hirogen_Functions.region_cutter(
                    shift, 
                    wavelengths,
                    flux*u.Unit('erg cm-2 s-1 AA-1'),
                    line_location - 40,
                    line_location + 40,
                    mode='Wavelength'
                )
            continuum = Hirogen_Functions.continua_maker(
                    wave_region,
                    flux_region,
                    # shift=Shift_Region,
                    item,
                    line_location,
                    object_name
                )
            peak_start = find_nearest(wavelengths, continuum[2][0])
            flux_without_peaks[peak_start:peak_start+len(continuum[0])] = continuum[0]

    return flux_without_peaks
    

def fits_file_gen(wave, flux, flux_err, flags, Main_Spectra_Path, Plate, MJD, FiberID, Survey, Run2D,
            Path_Override_Flag, Path_Override, Mode):

    hdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='flux', format='E', array = flux),
        fits.Column(name='loglam', format='E', array = wave),
        fits.Column(name = 'ivar', format = 'E', array = flux_err),
        fits.Column(name = 'and_mask', format = 'E', array = flags)]
    )

    if Mode.lower() == 'peaks':
        if Path_Override_Flag == 0:

            if Survey in ['sdss', 'segue1', 'segue2']:             

                filepath = f"{Main_Spectra_Path}/dr16/sdss/spectro/redux/{Run2D}/spectra/{Plate}/spec-{Plate}-{MJD}-{FiberID}-peaks.fits"
                                

            elif Survey in ['eboss', 'boss']:
            
                filepath = f"{Main_Spectra_Path}/dr16/eboss/spectro/redux/{Run2D}/spectra/full/{Plate}/spec-{Plate}-{MJD}-{FiberID}_peaks.fits"

            else:
                print("Not sure how to handle this survey - exiting for now")
                print(f"{Survey}")
                sys.exit()
        
        else:
            filepath = f"{Path_Override}"
    
    elif Mode.lower() == 'fakes':
        if Path_Override_Flag == 0:

            if Survey in ['sdss', 'segue1', 'segue2']:             

                filepath = f"{Main_Spectra_Path}/dr16/sdss/spectro/redux/{Run2D}/spectra/{Plate}/spec-{Plate}-{MJD}-{FiberID}-fake.fits"
                                

            elif Survey in ['eboss', 'boss']:
            
                filepath = f"{Main_Spectra_Path}/dr16/eboss/spectro/redux/{Run2D}/spectra/full/{Plate}/spec-{Plate}-{MJD}-{FiberID}_fake.fits"

            else:
                print("Not sure how to handle this survey - exiting for now")
                print(f"{Survey}")
                sys.exit()
        
        else:
            filepath = f"{Path_Override}"

    hdu.writeto(f"{filepath}", overwrite = True)
    
    return None

def find_nearest(array, number):
    array = np.asarray(array)
    if type(number) == u.quantity.Quantity:
        number = number.value
    idx = (np.abs(array - number)).argmin()
    return idx

def create_bins(sorted_pEQWs, bin_size):
    
    start = sorted_pEQWs[0]
    lower_lim = np.floor(start*2)/2

    end = sorted_pEQWs[-1]
    upper_lim = np.ceil(end*2)/2

    number_of_bins = int((np.abs(lower_lim) + np.abs(upper_lim))/bin_size)+1
    
    bins = np.linspace(lower_lim, upper_lim, number_of_bins, endpoint = True)
    return bins

def spectra_resizer(ecle_wavelengths, ecle_flux, galaxy_wavelengths, galaxy_flux):

    if ecle_wavelengths[0] > galaxy_wavelengths[0]:
        start_idx = find_nearest(galaxy_wavelengths, ecle_wavelengths[0])
        galaxy_wavelengths = np.delete(galaxy_wavelengths, np.arange(start_idx))
        galaxy_flux = np.delete(galaxy_flux, np.arange(start_idx))

    elif ecle_wavelengths[0] < galaxy_wavelengths[0]:
        start_idx = find_nearest(ecle_wavelengths, galaxy_wavelengths[0])
        ecle_wavelengths = np.delete(ecle_wavelengths, np.arange(start_idx))
        ecle_flux = np.delete(ecle_flux, np.arange(start_idx))

    if ecle_wavelengths[-1] < galaxy_wavelengths[-1]:
        end_idx = find_nearest(galaxy_wavelengths, ecle_wavelengths[-1]) + 1
        galaxy_wavelengths = np.delete(galaxy_wavelengths, np.arange(end_idx, len(galaxy_wavelengths)))
        galaxy_flux = np.delete(galaxy_flux, np.arange(end_idx, len(galaxy_flux)))

    if ecle_wavelengths[-1] > galaxy_wavelengths[-1]:
        end_idx = find_nearest(ecle_wavelengths, galaxy_wavelengths[-1]) + 1
        ecle_wavelengths = np.delete(ecle_wavelengths, np.arange(end_idx, len(ecle_wavelengths)))
        ecle_flux = np.delete(ecle_flux, np.arange(end_idx, len(ecle_flux)))

    return ecle_wavelengths, ecle_flux, galaxy_wavelengths, galaxy_flux