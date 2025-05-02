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

# Read in  CHARIS Data
filepath = '/home/ljd4475/projects/51erib/combined_data/data20170831/cubed_dataset/*.fits'
image_files = glob.glob(filepath)
image_files.extend(glob.glob('/home/ljd4475/projects/51erib/combined_data/data20170907/cubed_dataset/*.fits'))
image_files.extend(glob.glob('/home/ljd4475/projects/51erib/combined_data/data20180207/cubed_dataset/*.fits'))
image_files.extend(glob.glob('/home/ljd4475/projects/51erib/combined_data/data20190113/cubed_dataset/*.fits'))
image_files.sort()
dataset = CHARIS.CHARISData(image_files)

# list of null points 
points = [(90.05389998861241, 76.02299999938118),
          (109.92134452853566, 76.01274594264083),
          (123.97700000061882, 90.05389998861241),
          (123.98725405735917, 109.92134452853566),
          (109.94610001138759, 123.9770000006188),
          (90.07865547146434, 123.98725405735917),
          (76.0229999993812, 109.94610001138761),
          (76.01274594264083, 90.07865547146434)]

outpaths = [None,
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt1',
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt2',
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt3',
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt4',
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt5',
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt6',
            '/home/ljd4475/projects/51erib/MANUAL_nullpts/nullpt7']

prefixes = [None,
            'np1',
            'np2',
            'np3',
            'np4',
            'np5',
            'np6',
            'np7']

# setup FM guesses
# Used estimated position of x,y=(90,76) from interpolation of de Rosa paper.
# pyKLIP requests separation and angle so (x,y) were 
# converted to polar coordinates with the origin in the center (100,100) 
# and 0 degrees being upwards.
x = points[7][0]
y = points[7][1]
x_trans = x-100;
y_trans = y-100;
guesssep = np.sqrt(x_trans**2+y_trans**2)
guesspa = np.rad2deg(np.arctan2(y_trans,x_trans))-90
while guesspa < 0:
        guesspa += 360

guessflux = 1.49e+03 # In units of contrast to the host star 
                     # (put mouse over region)
star_type = 'F0V' # Found on Wikipedia

# Set up fm_class
wvs = np.unique(dataset.wvs) # in microns
dataset_gridamp = float(dataset.prihdrs[0]['X_GRDAMP'])
star_to_spot_ratio = 1./((dataset_gridamp 
                              / dataset.ref_spot_contrast_gridamp) ** 2
                              * dataset.ref_spot_contrast 
                              * dataset.ref_spot_contrast_wv ** 2 / (wvs ** 2))

# Spectral Extraction
exspec_outpath = outpaths[7]
prefix = prefixes[7] # fileprefix for the output files

# use the known planet separation and position angle,
# for example, use the measurements from the forward-model fitted astrometry
planet_sep = guesssep # companion separation in pixels
planet_pa = guesspa # companion position angle in degrees
planet_stamp_size = 15 # how big of a stamp around the companion in pixels, 
                       # stamp will be stamp_size**2 pixels
stellar_template = None # a stellar template spectrum, if you want

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

# Spectral Calibration
    
# obtain number of files and frames
N_frames = len(dataset.input)
N_cubes = np.size(np.unique(dataset.filenums))
nl = N_frames // N_cubes

# 3 spectral libraries available: 'ck04models', 'k93models', 'phoenix'
model_lib = 'ck04models'
temperature = 7331
metallicity = -0.12
log_g = 3.95
sp = pysynphot.Icat(model_lib, temperature, metallicity, log_g)
sp.convert("micron")
sp.convert("flam")
R_star = 1008765 # radius in km of 51 Eri A 
D_star = 9.082e+14 # distance in km to 51 Eri A - Wikipedia
spot_to_star_ratio = 1/star_to_spot_ratio
star_spectrum = np.array([sp.sample(dataset.wvs[i]) 
                          for i in range(nl)])*(R_star / D_star)**2
exspect_flux = exspect * star_spectrum * spot_to_star_ratio

# Save Spectra
i = 0
for spectrum in exspect_flux:
    wvs = list(wvs) 
    F = list(spectrum) 
    dict = {f'Wavelength ({sp.waveunits})': wvs, 
            f'F ({sp.fluxunits})': F}
    df = pd.DataFrame(dict)
    df.to_csv(exspec_outpath+f'/{prefix}_spectra_KL{numbasis[i]}.csv')
    i = i + 1
