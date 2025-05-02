# This script outputs PSF subtracted images of the data and extracts the 
# spectra for given KL values at a specified position and saves the spectra in 
# csv files. Before using this code, be sure pysynphot is installed properly.

# Import dependencies
import numpy as np
import glob
import pyklip.fm as fm
import pyklip.fmlib.extractSpec as es
from pyklip.instruments import CHARIS
import pysynphot
import pandas as pd
from uncertainties import ufloat
from scipy import stats as spstat
import pyklip.fakes as fakes

# Read in  CHARIS Data
filepath = '/home/ljd4475/projects/51erib/combined_data/data20170831/cubed_dataset/*.fits'
image_files = glob.glob(filepath)
image_files.sort()
dataset = CHARIS.CHARISData(image_files)

# setup FM guesses
# Used estimated position of x,y=(90,76) from interpolation of de Rosa paper.
# pyKLIP requests separation and angle so (x,y) were 
# converted to polar coordinates with the origin in the center (100,100) 
# and 0 degrees being upwards.
x = 90.0539
y = 76.023
x_trans = x-100;
y_trans = y-100;
guesssep = np.sqrt(x_trans**2+y_trans**2)
guesspa = np.rad2deg(np.arctan2(y_trans,x_trans))-90
while guesspa < 0:
        guesspa += 360
star_type = 'F0V' # Found on Wikipedia

# Set up fm_class
wvs = np.unique(dataset.wvs) # in microns
dataset_gridamp = float(dataset.prihdrs[0]['X_GRDAMP'])
star_to_spot_ratio = 1./((dataset_gridamp 
                              / dataset.ref_spot_contrast_gridamp) ** 2
                              * dataset.ref_spot_contrast 
                              * dataset.ref_spot_contrast_wv ** 2 / (wvs ** 2))

# Spectral Extraction
exspec_outpath = '/home/ljd4475/projects/51erib/all_4_dates/date1'
prefix = 'date1' # fileprefix for the output files

# use the known planet separation and position angle,
# for example, use the measurements from the forward-model fitted astrometry
planet_sep = guesssep # companion separation in pixels
planet_pa = guesspa # companion position angle in degrees
planet_stamp_size = 15 # how big of a stamp around the companion in pixels, 
                       # stamp will be stamp_size**2 pixels
stellar_template = 'methane' # a stellar template spectrum, if you want. Put None if not.

# reduction parameters
numbasis = np.array([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
maxnumbasis = 150
mode = 'ADI+SDI'
annuli=[[planet_sep-planet_stamp_size, planet_sep+planet_stamp_size]]
phi_section_size = 2 * planet_stamp_size / planet_sep # radians
subsections=[[np.radians(planet_pa) - phi_section_size / 2.,
              np.radians(planet_pa) + phi_section_size / 2.]]
movement = 2

# generate a background-subtracted satellite spot PSFs with shape 
# (nwv, boxsize, boxsize), averaged over exposures
boxsize = 15
dataset.generate_psfs(boxrad=boxsize//2)

fm_class = es.ExtractSpec(dataset.input.shape,
                          numbasis,
                          planet_sep,
                          planet_pa,
                          dataset.psfs,
                          np.unique(dataset.wvs),
                          stamp_size = planet_stamp_size)
fm.klip_dataset(dataset, fm_class,
                mode=mode,
                fileprefix=prefix,
                annuli=annuli,
                subsections=subsections,
                movement=movement,
                numbasis=numbasis,
                maxnumbasis=maxnumbasis,
                spectrum=stellar_template,
                outputdir=exspec_outpath, time_collapse='weighted-mean')

fmout_nanzero = np.copy(dataset.fmout)
fmout_nanzero[np.isnan(fmout_nanzero)] = 0.
exspect, fm_matrix = es.invert_spect_fmodel(fmout_nanzero, dataset,
                                            method='leastsq')


# Error Bar Calculation
def fake_spect(pa, fake_flux, basis):
    psfs = np.tile(dataset.psfs, (N_cubes, 1, 1))
    fakepsf = psfs * fake_flux[0,None,None]
    
    #tempdataset = charis.CHARISData(image_files,skipslices=[6,14,15,21])
    tempdataset = CHARIS.CHARISData(image_files)
    fakes.inject_planet(tempdataset.input, tempdataset.centers, fakepsf, 
                        tempdataset.wcs, planet_sep, pa)
    
    fm_class = es.ExtractSpec(tempdataset.input.shape,
                       basis,
                       planet_sep,
                       pa,
                       dataset.psfs,
                       np.unique(dataset.wvs),
                       stamp_size = planet_stamp_size)
    
    annuli=[[planet_sep-planet_stamp_size, planet_sep+planet_stamp_size]]
    phi_section_size = 2 * planet_stamp_size / planet_sep # radians
    subsections=[[np.radians(pa) - phi_section_size / 2.,
                  np.radians(pa) + phi_section_size / 2.]]

    fm.klip_dataset(tempdataset, fm_class,
                fileprefix="fmspect"+"pa_"+str(pa),
                annuli=annuli,
                subsections=subsections,
                movement=movement,
                numbasis = basis,
                spectrum=stellar_template,
                save_klipped=True, highpass=True,
                outputdir=exspec_outpath)
    
    exspect_fake, fm_matrix_fake = es.invert_spect_fmodel(tempdataset.fmout, 
                                            tempdataset,
                                            scaling_factor=None,
                                            method="leastsq")
    
    del tempdataset
    return exspect_fake

# obtain number of files and frames
N_frames = len(dataset.input)
N_cubes = np.size(np.unique(dataset.filenums))
nl = N_frames // N_cubes

npas = 11
pas = (np.linspace(planet_pa, planet_pa+360, num=npas+2)%360)[1:-1]
fake_spectra_all_bases = []
for i in range(len(numbasis)):
    input_spect = exspect[i,:]
    fake_spectra = np.zeros((npas, nl))
    for p, pa in enumerate(pas):
        fake_spectra[p,:] = fake_spect(pa, input_spect, numbasis[i])
    
    fake_spectra_all_bases.append(fake_spectra)

# extract error values
error = []
for i in range(len(fake_spectra_all_bases)):
    for j in range(nl):
        x = fake_spectra_all_bases[i][:,j]
        err = spstat.iqr(x)
        error.append(err)
e = np.array(error).reshape(int(len(error)/nl),nl)

# Spectral Calibration
    
# obtain number of files and frames
N_frames = len(dataset.input)
N_cubes = np.size(np.unique(dataset.filenums))
nl = N_frames // N_cubes

# 3 spectral libraries available: 'ck04models', 'k93models', 'phoenix'
model_lib = 'ck04models'
temperature = 7331
metallicity = -0.027
log_g = 3.95
sp = pysynphot.Icat(model_lib, temperature, metallicity, log_g)
sp.convert("micron")
sp.convert("flam")
# ufloat syntax: 1 plus/minus 2 = ufloat(1, 2)
R_star = ufloat(1008765, 13914) # radius in km of 51 Eri A - Wikipedia
# distance in km to 51 Eri A - Wikipedia
D_star = ufloat(9.0718920895846e+14,9257032744474.1) 
spot_to_star_ratio = 1/star_to_spot_ratio
star_spectrum = np.array([sp.sample(dataset.wvs[i]) 
                          for i in range(nl)])*(R_star / D_star)**2

#convert to ufloat
new_exspect = np.ndarray(shape=np.shape(exspect), dtype=object)
for i in range(len(exspect)):
    for j in range(nl):
        new_exspect[i][j] = ufloat(exspect[i,j],e[i,j])
        
#perform actual calibration
exspect_flux = new_exspect * star_spectrum * spot_to_star_ratio

#separate values
actual_values = np.zeros(np.shape(exspect))
error_bars = np.zeros(np.shape(exspect))
for i in range(len(exspect)):
    for j in range(nl):
        actual_values[i][j] = exspect_flux[i][j].nominal_value
        error_bars[i][j] = exspect_flux[i][j].std_dev

# Save Spectra
for i in range(len(numbasis)):
    dict = {f'Wavelength ({sp.waveunits})': dataset.wvs[:nl], # rather than wvs
            f'F ({sp.fluxunits})': actual_values[i],
            f'Error ({sp.fluxunits})': error_bars[i]}
    df = pd.DataFrame(dict)
    df.to_csv(exspec_outpath+f'/{prefix}_spec_KL{numbasis[i]}.csv')