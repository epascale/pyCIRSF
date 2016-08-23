import numpy as np
from scipy.interpolate import interp2d
import medfilt




def medfilt2d(ima, mask=None, wy = 3, wx = 3, step = 1, kind='cubic'):
    if not mask: mask = np.zeros(ima.shape, dtype=np.bool)
    
    if step <= 0: raise ValueError('invalid step value')
    
    ima_med = medfilt.medfilt2d(ima, mask, wy, wx, step)
    
    if step == 1: return ima_med
    
    # Interpolate over gaps using optimal interpolation
    idx0 = np.arange(0, ima.shape[0])
    idx1 = np.arange(0, ima.shape[1])
    
    print ima_med.shape, idx0.shape, idx0[::step].shape
    fp = interp2d(idx0[::step], idx1[::step], ima_med[::step,::step], kind=kind)
    
    return fp(idx0, idx1)

    


    

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
