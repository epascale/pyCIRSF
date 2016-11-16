import numpy as np
from astropy.io import ascii, fits
from astropy.table import hstack, Table, Column
import astropy.coordinates as coord
from astropy.stats import sigma_clip
import astropy.units as u
import scipy.interpolate
from photutils import CircularAperture, aperture_photometry, \
    CircularAnnulus, RectangularAnnulus
from openpyxl import load_workbook
import sys, os


def pycirsf_error(error_msg):
    sys.stderr.write("Error code: {:s}\n".format(error_msg))
    sys.exit(0)

def get_frames_flags(fname, date, object_name):
    wb = load_workbook(fname)
    
    ws = wb[date]
    
    ws_frames    = np.array([c.value for c in ws.columns[ 0]][1:])
    ws_objects   =          [c.value for c in ws.columns[ 1]][1:]
    ws_itimes    = np.array([c.value for c in ws.columns[ 2]][1:])
    ws_ra_off    = np.array([c.value for c in ws.columns[ 3]][1:])
    ws_de_off    = np.array([c.value for c in ws.columns[ 4]][1:])
    ws_flags     = np.array([c.value for c in ws.columns[14]][1:])
    ws_flags_new = np.array([c.value for c in ws.columns[15]][1:])
    
    list_ = []
    
    for i in xrange(len(ws_objects)):
	flag = np.int(ws_flags[i]) | np.int(ws_flags_new[i])
        if ws_objects[i] == object_name and flag == 0:
            list_.append([np.int(ws_frames[i]), ws_itimes[i], ws_ra_off[i], ws_de_off[i]])
        
    if list_:
        return Table([l for l in zip(*list_)], names=('Frame', 'ITIME', 'RA_OFF', 'DEC_OFF'))
    
    return list_
        
def get_reference_cat(fname=None, Jmag_lim=10.6):
    '''
    Read catalogue of star locations generated by starlink, 
    and retain stars brighter than J_maglim
    
    '''
    if not fname: pycirsf_error('Catalog file name not defined')
        
    cat = ascii.read(fname, header_start=14, data_start=16)
    cat.sort('Jmag')
    
    idx = np.where(cat['Jmag'] < Jmag_lim)
    cat = cat[idx]
    
    return cat, coord.Angle(cat['RAJ2000'], u.hour), \
        coord.Angle(cat['DEJ2000'], u.degree)

def get_dark(fname, flag_mask = 0x1|0x4):
    '''
    Read dark frame from fname and returns a masked array
    '''
    
    if not os.path.exists(fname): pycirsf_error('Mask file ({:s}) does not exist'.format(fname))
    hdu = fits.open(fname)
    ima  = hdu[1].data
    mask = hdu[2].data
    hdu.close()
    
    mask = np.bitwise_and(mask, flag_mask).astype(np.bool)
    
    return np.ma.array(ima, mask=mask)

def get_flat(fname, flag_mask = 0x1|0x4):
    '''
    Read dark frame from fname and returns a masked array
    '''
    
    if not os.path.exists(fname): pycirsf_error('Mask file ({:s}) does not exist'.format(fname))
    hdu = fits.open(fname)
    ima  = hdu[0].data
    
    ima = sigma_clip(ima, sigma = 5.0)
    
    #mask = hdu[2].data
    hdu.close()
    
    #mask = np.bitwise_and(mask, flag_mask).astype(np.bool)
    
    return ima
    
def apply_dark_flat(ima, dark=None, flat=None, fillval = 0.0):
    
    ima_ = ima
    
    if isinstance(dark, np.ndarray) : ima_ = ima_ - dark
    if isinstance(flat, np.ndarray):  ima_ = ima_ / flat
    
    ima_.set_fill_value(fillval)
    
    return ima_
 
def multiplicative_mask(mask, flag_mask = 0x1|0x4, nanfill=None):
    ''' 
    Create a mask array containing 0 (bad pixels) and 1 (good pixels)
    from the raw mask array. 
    Mask values are specidied in the mask header and currently are
    0x1: Array central region
    0x2: |QE-1| > 3sigma
    0x4: |QE-1| > 5sigma
    0x8: |QE-1| > 10sigma
    
    '''
    if not np.issubdtype(mask.dtype, np.integer):
        pycirsf_error('multiplicative_mask: mask dtype should be int')
    
    mask_ = np.bitwise_xor(np.bitwise_and(mask, flag_mask), flag_mask)/flag_mask
    
    if nanfill:
      idx = np.where(mask_ == 0)
      mask_ = mask_.astype(np.float32)
      mask_[idx] = np.nan
      
    return mask_
    
def masked_array(ima, mask, flag_mask = 0x1|0x4):
    ''' 
    Create a masked array from the raw image and raw mask array.
    See pyCRSF.lib.multiplicative_mask for availabe mask flags, or mask fits header

    '''

    if not np.issubdtype(mask.dtype, np.integer):
        pycirsf_error('multiplicative_mask: mask dtype should be int')
    
    mask_ = np.bitwise_and(mask, flag_mask).astype(np.bool)
    ima = ima*np.logical_not(mask_).astype(np.float)
    return np.ma.masked_array(ima, mask_)

    

def remove_bias(ima, mask, flag_mask = 0x1|0x4, region='quadrant', zero_off=False):
    '''
    Remove bias from image by estimating the median of the array 
    the region-of-interest specified by 'region'.
    
    Parameters
    ----------
    
    ima  : array_like
           input image to be de-biased
    mask : array_like
           mask to be used
    flag : scalar
           flag value to apply in mask. Mask pixels that match flag in the bitwise or are excluded.
    region : string
           if 'quadrant' is specified, the array is divided in 4 quadrants
           and a median value per quadrant is estimated.
           If 'section' is specified, the array is divided into two sections (left and write
           along axis 0) and a median is estimated in each of the two sections.
    zero_off: bool
           if False (default) each section is scaled such that its DC-level metches the average 
           of medians estimated. If True, the median in each quadrant is subtracted from each quadrant.
    
    '''
    
    dy, dx = ima.shape
    
    maskm = multiplicative_mask(mask, flag_mask=flag_mask)
    
    m = []
    if region == 'quadrant':
        pass
    
    elif region == 'section':
        i = 1
        for j in [0, 1]:
           idx = np.where(maskm[..., j*dx/2:(j+1)*dx/2])
           m.append(np.median(ima[..., j*dx/2:(j+1)*dx/2][idx]))
        
        for j in [0, 1]:
            ima[..., j*dx/2:(j+1)*dx/2] = ima[..., j*dx/2:(j+1)*dx/2] - m[i]
    
    return ima
           


def stacking_(ima, xc, yc, N = 31, remove_background=True):
    '''
    Stack the input image at locations xc, yc. 
    The stack image has shape NxN. 
    If remove_background is true, a background estimated in a rectangular annulus and removed.
    
    '''
    
    ima_ = ima.copy()
    idx = np.where(np.isnan(ima_))
    ima_[idx] = 0.0

    stamp = np.zeros( (N, N) )
    
    x = np.arange(ima_.shape[1])
    y = np.arange(ima_.shape[0])
    
    f = scipy.interpolate.interp2d(x, y, ima_, kind='cubic', fill_value=0.0)
    
    new_x = np.arange(N) - N//2
    new_y = np.arange(N) - N//2
    
    for xx, yy in zip(xc, yc):
        stamp += f(new_x+xx, new_y+yy)
    
    if remove_background:
        anulus_apertures = RectangularAnnulus(
            (stamp.shape[1]/2, stamp.shape[0]/2), 
            w_out=stamp.shape[0], h_out=stamp.shape[1],
            w_in=0.8*stamp.shape[0], theta=0.0)
        
        bkg = aperture_photometry(stamp, anulus_apertures)
        bkg['aperture_sum'] /= anulus_apertures.area()
        stamp -= bkg['aperture_sum']
        
    
    stamp /= stamp.sum()
    
    light_baricenter_x = (stamp*np.arange(stamp.shape[1])).sum()
    light_baricenter_y = (stamp.transpose()*np.arange(stamp.shape[0])).sum()
    
    return stamp, light_baricenter_x, light_baricenter_y

def stacking(ima, xc, yc, N = 31, remove_background=True, method='mean'):
    '''
    Stack the input image at locations xc, yc. 
    The stack image has shape NxN. 
    If remove_background is true, a background estimated in a rectangular annulus and removed.
    
    '''
    
    ima_ = ima.copy()

    stamp = []
    for xx, yy in zip(xc, yc):
        xx_ = np.int(xx); yy_ = np.int(yy)
        stamp_ = ima_[-N//2+1+yy_:N//2+1+yy_, -N//2+1+xx_:N//2+1+xx_]
        if stamp_.shape[0] == N and stamp_.shape[1] == N :
	  stamp.append( stamp_ )
    
    if method=='mean':
      stamp = np.ma.array(stamp).mean(axis=0)
    elif method=='median':
      stamp = np.ma.median(np.ma.array(stamp), axis=0)
    else:
      pycirsf_error('Stacking: method not valid')
      
    if remove_background:
        anulus_apertures = RectangularAnnulus(
            (stamp.shape[1]/2, stamp.shape[0]/2), 
            w_out=stamp.shape[0], h_out=stamp.shape[1],
            w_in=0.8*stamp.shape[0], theta=0.0)
        
        bkg = aperture_photometry(stamp, anulus_apertures)
        bkg['aperture_sum'] /= anulus_apertures.area()
        stamp -= bkg['aperture_sum']
        
    
    stamp /= stamp.sum()
    
    light_baricenter_x = (stamp*np.arange(stamp.shape[1])).sum()
    light_baricenter_y = (stamp.transpose()*np.arange(stamp.shape[0])).sum()
    sx = stamp.sum(axis=0)
    sx /= sx.sum()
    light_baricenter_x = (sx*np.arange(stamp.shape[1])).sum()
    sy = stamp.sum(axis=1)
    sy /= sy.sum()
    light_baricenter_y = (sy*np.arange(stamp.shape[0])).sum()
    
    return stamp, light_baricenter_x, light_baricenter_y

def stacking_nolb(ima, xc, yc, N = 31, remove_background=False):
    '''
    Stack the input image at locations xc, yc. 
    The stack image has shape NxN. 
    If remove_background is true, a background estimated in a rectangular annulus and removed.
    
    '''
    
    ima_ = ima.copy()

    # Need to ensure that sources are not too close to edge
    badx = [np.where((xc<N/2.) | (xc>np.shape(ima)[0] - N/2.))]
    xc = np.delete(xc, badx)
    yc = np.delete(yc, badx)
    bady = [np.where((yc<N/2.) | (yc>np.shape(ima)[0] - N/2.))]
    xc = np.delete(xc, bady)
    yc = np.delete(yc, bady)

    stamp = []
    for xx, yy in zip(xc, yc):
        xx_ = np.int(xx)
        yy_ = np.int(yy)
        stamp.append( ima_[-N//2+1+yy_:N//2+1+yy_, -N//2+1+xx_:N//2+1+xx_] )
    
    stamp = np.ma.array(stamp).mean(axis=0)
    
    if remove_background:
        anulus_apertures = RectangularAnnulus(
            (stamp.shape[1]/2, stamp.shape[0]/2), 
            w_out=stamp.shape[0], h_out=stamp.shape[1],
            w_in=0.8*stamp.shape[0], theta=0.0)
        
        bkg = aperture_photometry(stamp, anulus_apertures)
        bkg['aperture_sum'] /= anulus_apertures.area()
        stamp -= bkg['aperture_sum']
        
    
    stamp /= stamp.sum()
    
    return stamp

def barycenter(image):
    """
    This function will calculate the lightbarycenter and the expected values of
    x**2, y**2, x**3 and y**3 for a given input image (or stamp).
    """

    # Pixel values covering range of stamp
    xx = np.arange(image.shape[1])
    yy = np.arange(image.shape[0])
    
    # Need y-coords of every pixel 
    pointsy = np.zeros(len(xx)**2)
    for i in xx:
        g = i*len(xx)
        pointsy[g:g+len(xx)] = i
    
    # Need x-coords of every pixel (needs to correspond to order of y)
    pointsx = np.zeros(len(yy)**2)
    for j in yy:
        h = j*len(xx)
        pointsx[h:h+len(xx)] = yy
    
    # Need to flatten image   
    data = image.flatten()
    
    # Light barycenter       
    lbx = ((data*pointsx).sum()) / (data.sum())
    lby = ((data*pointsy).sum()) / (data.sum())
    
    # Expected value for variance
    xsqu = ((data*(pointsx**2)).sum()) / (data.sum())
    ysqu = ((data*(pointsy**2)).sum()) / (data.sum())
    
    # Expected value for asymmetry
    xcub = ((data*(pointsx**3)).sum()) / (data.sum())
    ycub = ((data*(pointsy**3)).sum()) / (data.sum())     
  
    return lbx, lby, xsqu, ysqu, xcub, ycub 

def moments(image):
  
  norm = 1.0/image.sum()
  xc = np.sum(image * np.arange(image.shape[1]))*norm
  yc = np.sum(image.transpose() * np.arange(image.shape[0]))*norm
     
  return yc, xc

def photom(ima, pos, radius, r_in=False, r_out=False, mode='median'):
    '''
    Aperture photometry in an aperture located at pixel coordinates 
    pos = ( (x0, y0), (x1, y1), ... ) with a radius=radius.
    When r_in and r_out are given, background is estimated in CircularAnnulus and subtracted.
    
    mode refers to how the background is estimated within the circlar annulus.
    Can be 'median' or 'mean'

    '''
    # Setting up the mask 
    if hasattr(ima, 'mask'):
      if ima.mask.size == 1:
	mask = np.zeros(ima.shape, dtype=np.bool) | ima.mask
      else:
        mask = ima.mask.copy()
    else:
        mask = np.zeros(ima.shape, dtype=np.bool)
    ### Performing the actual photometry - identical for each method
    # Setting up the aperture 
    apertures = CircularAperture(pos, r = radius)
    # Aperture photometry on image
    ap        = aperture_photometry(ima, apertures, mask=mask)
    # Aperture photometry on mask to see how many masked pixels are in the 
    # aperture
    apm       = aperture_photometry(mask.astype(int), apertures)
    # Number of masked pixels in aperture
    map_area  = Column(name='bpix_aper', data= apm['aperture_sum'].data)   
    # Number of unmasked pixels in aperture
    ap_area   = Column(name = 'area_aper',
			    data=apertures.area() - apm['aperture_sum'].data)
    # Flux of pixels
    flux_init      = Column(name = 'flux', data=ap['aperture_sum'].data)
        
    bkg = False
    ### Two different modes for analysing the background
    if ( r_in and r_out and mode in ('mean', 'median') ):
      ### This stuff is the same regardless of method
      # Setting up the annulus
      anulus_apertures = CircularAnnulus(pos, r_in=r_in, r_out=r_out)
      # Performing annulus photometry on the mask
      bkgm = aperture_photometry(mask.astype(int), anulus_apertures)
      # Number of masked pixels in bkg
      mbkg_area = Column(name = 'bpix_bkg',
			 data=bkgm['aperture_sum'])  
      # Number of non-masked pixels in aperture and bkg        
      bkg_area  = Column(name = 'area_bkg',
			 data=anulus_apertures.area() - bkgm['aperture_sum'])
      # Adding this data to table
      ap.add_column(bkg_area)
      ap.add_column(mbkg_area)
      
      ### This stuff is specific to the mean
      if mode == 'mean':
	# Perform the annulus photometry on the image
	bkg  = aperture_photometry(ima, anulus_apertures, mask=mask)
        # Average bkg where this divides by only number of NONMASKED pixels
        # as the aperture photometry ignores the masked pixels
        bkga = Column(name='background',
		      data=bkg['aperture_sum']/bkg_area)
        # Bkg subtracted flux
        flux = flux_init - bkga*ap_area
        # Adding that data
        ap.add_column(bkga)
      elif mode == 'median':
	# Number of pixels in the annulus, a different method
	fractions = anulus_apertures.get_fractions(ima, method='center')
	nbkg = fractions.shape[-1] if fractions.ndim == 3 else 1
	# Background mask
	bkgm = np.zeros(nbkg, dtype=np.float)
	
	if bkgm.size == 1:
	  bmask = ~mask & fractions.astype(np.bool)
	  bkgm[0] = np.median(ima[bmask])
	else:	
	  for i in xrange(bkgm.size):
	    bmask = ~mask & fractions[..., i].astype(np.bool)
	    bkgm[i] = np.median(ima[bmask])
		
	flux = flux_init - bkgm*ap_area
	bkgm = Column(name = 'background', data = bkgm)
	ap.add_column(bkgm)
          
    ap.add_column(ap_area)
    ap.add_column(map_area)
    ap.add_column(flux)
        
    return ap['flux'], ap, map_area

# Find the nearest value to a number, return x and y 
def find_nearest(ydata, xdata, yvalue):
    # Trunkate ydata so that sigma value only considered before peak
    trunk = np.argmax(ydata)
    ydat = ydata[0:trunk]
    idx = (np.abs(ydat-yvalue)).argmin()
    return idx


class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.02f}'.format(x, y, z)



if __name__ == "__main__":
    from astropy.io import fits
    import matplotlib.pyplot as plt
    
    
    mask = fits.open(os.path.expanduser('~/IRSF/calibration/kmask_160128.fits'))[0].data
    flat = fits.open(os.path.expanduser('~/IRSF/calibration/kflat_160128.fits'))[0].data
    #flat = fits.open(os.path.expanduser('~/IRSF/pyIRSF-2.2.1/lib/kcflat.fits'))[0].data
    ima  = fits.open(os.path.expanduser('~/IRSF/data/160128/rawdata/k160128_0100.fits.fz'))[1].data
    
    maskm = multiplicative_mask(mask, nanfill=False)
    ima = maskm * ima / flat 
    #ima = remove_bias(ima, mask, flag_mask = 0x1|0x8, region='section')
    
    plt.ion()
    fig, (ax0, ax1) = plt.subplots(ncols=2, num=1)
    im0 = ax0.imshow(maskm, interpolation='None')
    ax0.format_coord = Formatter(im0)
    vlo, vhi = np.percentile(ima, (5, 99.99))
    im1 = ax1.imshow(ima, vmin=vlo, vmax=vhi, interpolation='None')
    ax1.format_coord = Formatter(im0)
    
    pass