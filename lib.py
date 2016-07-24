import numpy as np
import sys, os

def pycirsf_error(error_msg):
    sys.stderr.write("Error code: {:s}\n".format(error_msg))
    sys.exit(0)

 
 
def multiplicative_mask(mask, flag_mask = 0x1|0x8, nanfill=None):

    if not np.issubdtype(mask.dtype, np.integer):
        pycirsf_error('multiplicative_mask: mask dtype should be int')
    
    mask_ = np.bitwise_xor(np.bitwise_and(mask, flag_mask), flag_mask)/flag_mask
    
    if nanfill:
      idx = np.where(mask_ == 0)
      mask_ = mask_.astype(np.float32)
      mask_[idx] = np.nan
      
    return mask_
    

def remove_bias(ima, mask, flag_mask = 0x1|0x8, region='quadrant', zero_off=False):
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