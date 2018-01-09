import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel, Box2DKernel, Tophat2DKernel

from photutils import CircularAperture, aperture_photometry, CircularAnnulus
from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar

import pyCIRSF as irsf
import os, glob, sys, mylib


#######################################################
def photom(ima, pos, radius, r_in=None, r_out=None, method='median'):
    '''
    Aperture photometry in an aperture located at pixel coordinates 
    pos = ( (x0, y0), (x1, y1), ... ) with a radius=radius.
    When r_in and r_out are given, background is estimated in CircularAnnulus and subtracted.
    
    method refers to how the background is estimated within the circlar annulus.
    Can be 'median' or 'mean' or 'mode'

    '''
    
    ima_local = np.ma.asanyarray(ima.copy())
    ima_local.fill_value = np.nan
    mask_ = ima_local.mask
    ima_  = ima_local.filled()
    
    ### Do photometry - identical for each method
    apertures = CircularAperture(pos, r = radius)
    ap        = aperture_photometry(ima_, apertures, 
                                    mask=mask_, method='exact')
    # Aperture photometry on mask to estimate # of masked pixels in aperture
    apm       = aperture_photometry(mask_.astype(int), apertures,
                                    method='exact')
    ap.add_columns( [apertures.area()-apm['aperture_sum'], apm['aperture_sum']],
                   names=['aperture_area', 'aperture_badpix'])
    ap.add_column(ap['aperture_sum'], index=3, name='Flux')

    if ( r_in == None or r_out == None or not method in ('mean', 'median', 'mode') ): 
      # Quit here if background correction is not requested
      return ap

    annulus_apertures = CircularAnnulus(pos, r_in=r_in, r_out=r_out)
    annulus_masks = annulus_apertures.to_mask(method='center')
    bg_values = []
    for annulus_mask in annulus_masks:
      bg_ima = annulus_mask.cutout(ima_)
      bg_mask = annulus_mask.cutout(mask_.astype(np.int))
      bg_ima = np.ma.array(bg_ima, mask= bg_mask.astype(np.bool) | ~annulus_mask.data.astype(np.bool))
      if method == 'mean': bg_val = bg_ima.mean()
      elif method == 'median': bg_val = np.ma.median(bg_ima)
      elif method == 'mode': 
        kernel = gaussian_kde(bg_ima.data[~bg_ima.mask], bw_method='silverman')
        mode = bg_ima.mean()
        std  = bg_ima.std()
        
        mode = minimize_scalar(lambda x: -kernel(x), bounds=(mode-3*std, mode+3*std),
                               method='bounded')
        bg_val=mode.x[0]
        
        if True:
          median = np.ma.median(bg_ima)
          h, b = np.histogram(bg_ima.data[~bg_ima.mask], bins=15, normed=True)
          bc = 0.5*(b[1:]+ b[:-1])
          plt.figure(33); plt.clf(); plt.ioff()
          fig, (ax0,ax1) = plt.subplots(ncols=2, nrows=1, num=33)
          ax0.plot(bc, h, 'x')
          x = np.linspace(bc.min(), bc.max(), 100)
          ax0.plot(x, kernel(x))
          ax0.vlines(mode.x, ax0.get_ylim()[0], ax0.get_ylim()[1])
          ax0.vlines(median, ax0.get_ylim()[0], ax0.get_ylim()[1])
          ax1.imshow(bg_ima)
          plt.show()
        
        
      bg_values.append(bg_val)
    ap.add_column(Column(data=bg_values, name = 'background'))  
    ap['Flux'] = ap['Flux'] - ap['aperture_area']*ap['background']
    return ap, bg_ima
##################################################################

dates       = ['170420', '170428'][:1]
band        = 'k'
flat_type   = 'c'
object_name = 'WASP-103'
file_open_mode = 'w'
file_format = 'ascii'

__root__ = "~/gdrive/WorkV2/IRSF/Observations 2017"
workbook_fn = os.path.join(__root__,  'obslog_2017_jeni.xlsx')

flat_fn = os.path.join('data/flats/new_flats', '{:s}{:s}flat.fits'.format(band, flat_type))
flat = irsf.lib.get_flat(os.path.join(__root__, flat_fn), sigma=5.0)

fignum = 31415; plt.ion(); fig = plt.figure(fignum); plt.clf()
fig, ax0 = plt.subplots(nrows = 1, ncols = 1, num=fignum)


for date in dates:
#1# read frames
  raw_data_path = os.path.join(__root__, 'data', 'raw', date, 'rawdata')
  fr_tab = irsf.lib.get_frames_flags(os.path.expanduser(workbook_fn), 
                                     date, object_name)
#2# read dark frame
  for which_dark in ('as', 'bs'):    
    dark_fn = os.path.join(__root__, 'data/darks', object_name,
                          date,
                          'dark_masked_{:s}_{:s}_{:s}.fits.fz'.format(band, date, which_dark))
    if os.path.isfile(os.path.expanduser(dark_fn)): 
      dark = irsf.lib.get_dark(dark_fn, 
                               flag_mask=0x1, dark_ext = 1, mask_ext = 2)
  
      break
    
  for fr in fr_tab['Frame'][:3]:  
    print "Frame {:04d}\r".format(fr),
#3# read sci frame  
    raw_fn = os.path.expanduser(os.path.join(raw_data_path,
              '{:s}{:s}_{:04d}.fits.fz'.format(band, date, fr)))
    
    hdulist = fits.open(raw_fn) 
    hdr = hdulist[1].header
    ima = irsf.lib.apply_dark_flat(np.ma.array(hdulist[1].data), 
                  dark=dark, flat=flat, linearity_correction=True, band=band)
    hdulist.close()

#4# interpolate over hot pixels
#    kernel = Box2DKernel(3, mode='center')
#    ima_ = interpolate_replace_nans(ima.filled(), kernel)
#    ima = np.ma.array(ima_, mask=np.isnan(ima_))
    #dsds
#4# read positions
    sourcelist_fn = os.path.expanduser(os.path.join(__root__, "data/sources",
                                 "{:s}{:s}_{:04d}_sources.dat".format(band, date, fr)))
    
    if not os.path.isfile(sourcelist_fn): continue
    source_list = ascii.read(sourcelist_fn, header_start=0, comment='=')
    
    positions = np.asanyarray(source_list['X_CENTER_Gauss', 'Y_CENTER_Gauss']).tolist()
    photometry, itmp = photom(ima, positions, 5.0, r_in = 10, r_out=20, method='mode')
    for source in source_list:
      with open(source['2MASS'], file_open_mode) as fs:
        keys = ['2MASS', 'X_CENTER_Gauss', 'Y_CENTER_Gauss']
        names = ['2MASS', 'X', 'Y']
        output_table = Table([[source[key]] for key in keys], names=names)
        output_table.write(fs, format=file_format, delimiter=',')

    file_open_mode = 'a'
    file_format = 'ascii.no_header'

# # combile table
    
    

    
    vmin, vmax = np.percentile(ima.flatten(), (5, 99.))
    im = ax0.imshow(ima.filled(), interpolation='none', cmap='gist_heat', 
                   vmin=vmin, vmax=vmax, origin='lower')
    ax0.format_coord = irsf.lib.Formatter(im)
    ax0.grid()

  
  
