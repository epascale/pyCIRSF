import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.wcs as wcs

from photutils.background import Background2D
from photutils import CircularAperture
from astropy.convolution import convolve, Gaussian2DKernel
from scipy.signal import medfilt2d, medfilt

import pyCIRSF as irsf

import os, glob, sys


date = '160127'
band = 'h'
flat_suffix = 'flat_160128_v2'
flat_suffix = 'cflat'
dark_suffix = 'dark_master'
object_name = 'WASP-121'

Nstack = 127
r0 = 15

astrom_data_path = os.path.join('~/IRSF/proc_data_ep', date)
raw_data_path = os.path.join('~/IRSF/data', date, 'rawdata')
cal_path = '~/IRSF/calibration'
cat_fname = '~/IRSF/etc/ref_cat_2mass.dat'
wb_fname         = '~/IRSF/etc/IRSF_newflags_v2.xlsx'

cat_2mass, cat_ra, cat_dec = irsf.lib.get_reference_cat(
    fname=os.path.expanduser(cat_fname),
    Jmag_lim=10.6)

dark = irsf.lib.get_dark(os.path.expanduser(os.path.join(
    cal_path, band+dark_suffix+'.fits.fz')), flag_mask=0x1 | 0x4)

flat = irsf.lib.get_flat(os.path.expanduser(os.path.join(
    cal_path, band+flat_suffix+'.fits')), flag_mask=0x1 | 0x4)

fr_tab = irsf.lib.get_frames_flags(os.path.expanduser(wb_fname), 
				   date, object_name)


#plt.ion(); 
plt.figure(123); plt.clf()
fig, ( (ax0, ax1, ax2), (ax3, ax4, ax5) ) = plt.subplots(nrows=2, 
							 ncols=3, 
							 num=123)

for fr in fr_tab['Frame'][:1]:

    wcs_fn = os.path.expanduser(os.path.join(astrom_data_path,
            '{:s}{:s}_{:04d}.fits'.format(band, date, fr)))
    raw_fn = os.path.expanduser(os.path.join(raw_data_path,
            '{:s}{:s}_{:04d}.fits.fz'.format(band, date, fr)))

    hdulist = fits.open(wcs_fn)
    hdr     = hdulist[0].header
    hdulist.close()
    
    hdulist = fits.open(raw_fn)
    
    ima =irsf.lib.apply_dark_flat(hdulist[1].data, dark=dark, flat=flat)
    
    hdulist.close()
    
    w = wcs.WCS(hdr)
    
    
    print 'filtering ...'
    cat_2mass_, ra_, dec_ = irsf.lib.get_reference_cat(
	fname=os.path.expanduser(cat_fname),
	Jmag_lim=18)
    cat_x_, cat_y_ = w.all_world2pix(ra_.degree, dec_, 0)         
    ima, background = irsf.medianfilter.remove_background(ima, source_x = cat_x_, source_y = cat_y_, source_r = r0)
    print 'done.'
    
    # Create the stack
    cat_x, cat_y = w.all_world2pix(cat_ra.degree, cat_dec, 0)
    stack, lbx, lby = irsf.lib.stacking(ima, cat_x, cat_y, N=Nstack, 
					remove_background=False)
    
    pos = (lbx, lby)
    radii = np.linspace(1, 50, 40)
    flux = []
    for r in radii:
        flux_, ap = irsf.lib.photom(stack, pos, r, r_in=50, r_out=60)
        flux.append(flux_)
	
    vmin, vmax = np.percentile(ima.filled().flatten(), (5, 99.9))
    im0 = ax0.imshow(ima, interpolation='none', cmap='gist_heat', vmin=vmin, vmax=vmax)
    ax0.contour(background, colors='w', alpha=0.2)
    ax0.format_coord = irsf.lib.Formatter(im0)
    ax0.autoscale(False)
    ax0.plot(cat_x, cat_y, '+r')
    apertures = CircularAperture( (cat_x, cat_y), r = r0)
    apertures.plot(ax = ax0, color='r')
    
    im1 = ax1.imshow(np.log10(stack-stack.min()+1e-6), interpolation='none', cmap='gist_heat')
    ax1.format_coord = irsf.lib.Formatter(im1)
    
    ax2.plot(np.ma.median(ima, axis=0), 'r', label='gradient across x')
    ax2.plot(np.ma.median(ima, axis=1), 'b', label='gradient across y')
    #ax2.plot(np.ma.median(background, axis=0), 'm')
    #ax2.plot(np.ma.median(background, axis=1), 'c')
    ax2.legend()
    
    ax3.plot(radii, flux, 'o-k')
    ax3.grid()
    ax3t = ax3.twinx()
    ax3t.plot(stack.sum(axis=0)[stack.shape[1]//2:], 'r', label='$\int\, PSF\, dy$')
    ax3t.plot(stack.sum(axis=1)[stack.shape[0]//2:], 'b', label='$\int\, PSF\, dx$')
    
    ax4.plot(stack.sum(axis=0), 'r', label='$\int\, PSF\, dy$')
    ax4.plot(stack.sum(axis=1), 'b', label='$\int\, PSF\, dx$')
    ymin, ymax = ax4.get_ylim()
    ax4.vlines(lbx, ymin, ymax, colors='r')
    ax4.vlines(lby, ymin, ymax, colors='b')    
    ax4t = ax4.twinx()
    ax4t.plot(stack.sum(axis=0).cumsum(), 'r')
    ax4t.plot(stack.sum(axis=1).cumsum(), 'b')
    ax4.grid()
    ax4.legend()
    
    vmin, vmax = np.percentile(background.flatten(), (5, 99.9))
    im5 = ax5.imshow(background.filled(), interpolation='none', cmap='gist_heat', vmin=1000, vmax=1100)
    ax5.format_coord = irsf.lib.Formatter(im5)
    ax5.grid()

plt.show()    
    