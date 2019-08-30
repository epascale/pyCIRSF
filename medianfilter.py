import numpy as np
from scipy.interpolate import interp2d, RectBivariateSpline
import pyCIRSF.medfilt as medfilt





def medfilt2d(ima, mask=None, wy = 3, wx = 3, step = 1, kind='cubic'):
    # Copy image and set up mask if none there
    
    if hasattr(ima, "filled"):
      _ima = ima.filled().copy()
    else:
      _ima = ima.copy()
    
    _mask = np.zeros(_ima.shape, dtype=np.bool)
    
    if hasattr(ima, 'mask'):
      _mask = np.logical_or(_mask, ima.mask)
    
    if not mask is None:
      _mask = np.logical_or(_mask, mask)
      
      
    if step <= 0: raise ValueError('invalid step value')
    # Perform median filtering (reduces noise)
    ima_med = medfilt.medfilt2d(_ima, _mask, wy, wx, step)
    if mask is None and step == 1: return ima_med
    
    # Interpolate over gaps using optimal interpolation
    idx0 = np.arange(0, _ima.shape[0])
    idx1 = np.arange(0, _ima.shape[1])
    
    #fp = RectBivariateSpline(idx1[::step], idx0[::step], ima_med[::step,::step])
    
    # Removing inf
    idx = np.where(np.isinf(ima_med))
    if idx: ima_med[idx] = 0.0
    
    # Interpolate over median filtering to model background
    fp = interp2d(idx1[::step], idx0[::step], ima_med[::step,::step], kind=kind)
    
    return fp(idx1, idx0)

def remove_background(ima, source_x=None, source_y=None, source_r=None, wx=64, wy=64, step=8):
    # Make a copy of image and check for mask
    _ima = ima.copy()

    if hasattr(_ima, 'mask'):
      if _ima.mask.size == 1:
        _mask = np.zeros(ima.shape, dtype=np.bool) | ima.mask
        _ima = np.ma.array(_ima, mask=_mask)
    else:
        _mask = np.zeros(ima.shape, dtype=np.bool)
        _ima = np.ma.array(_ima, mask=_mask)
    
    if not (source_x is None or source_y is None or source_r is None):
      source_mask = np.zeros( _ima.shape, dtype=bool )
      
      for x, y in zip(np.round(source_x).astype(int), np.round(source_y).astype(int)):
        source_mask[y - source_r:y+source_r+1,x - source_r:x+source_r+1] = True
        _ima.mask |= source_mask
    
    # Model background by splitting image into two halves (L and R)  
    nx = ima.shape[1]
    bkg1 = medfilt2d(_ima.filled()[..., :nx//2], mask=_ima.mask[..., :nx//2], 
                                        wx=wx, wy=wy, step=step)
    bkg2 = medfilt2d(_ima.filled()[..., nx//2:], mask=_ima.mask[..., nx//2:], 
                                        wx=wx, wy=wy, step=step)
    # Put halves back together
    background = np.ma.array(np.hstack( (bkg1, bkg2) ), 
                              mask=ima.mask,
                              fill_value=ima.fill_value) 
    return ima-background, background
  

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    import scipy.signal
    
    d = 2
    wy, wx = 32/d+1, 32/d+1
    step = 8/d
    N = 1024/d

    ima = np.random.randn(N,N)
    mask = np.zeros( (N, N), dtype=bool)
    ima_med = medfilt2d(ima, wy=wy, wx=wx, step=step)
    ima_med_sp = scipy.signal.medfilt2d(ima, kernel_size=wx)


    
    
    #plt.ion(); 
    plt.figure(4); plt.clf()
    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, num=4)

    ax0.imshow(ima, interpolation='none', vmin=ima.min(), vmax=ima.max())
    ax1.imshow(ima_med, interpolation='none', vmin=ima_med.min(), vmax=ima_med.max())
    ax3.imshow(ima_med_sp, interpolation='none', vmin=ima_med_sp.min(), vmax=ima_med_sp.max())
    
    ima_diff = ima_med_sp-ima_med
    vmin, vmax = np.percentile(ima_diff, (1.0, 99.0))
    ax2.imshow(ima_diff, interpolation='none', vmin=vmin, vmax=vmax)
    plt.show()
