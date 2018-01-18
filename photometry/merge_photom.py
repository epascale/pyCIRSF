from astropy.io import ascii
import numpy as np
from astropy.table import Table, Column
import os

#######################################################
file_format = 'csv'
band = 'k'

r_ap = 3.5
method = 'mode'

r_ap_list = [3.5,4.5, 5.0, 10]
method_list = ['mode', 'median']


__root__ = "~/Documenti/gdrive/WorkV2/IRSF/Observations 2017/"

photom_dir = "data/photometry_180115"

sourcelist_fn = os.path.expanduser(os.path.join(__root__, "data/sources",
                             "{:s}{:s}_{:04d}_sources.dat".format(band, '170502', 381)))

source_list = ascii.read(sourcelist_fn, header_start=0, comment='=')

ids = source_list['2MASS']

for method in method_list:
    print '\nMethod ' + method

    for r_ap in r_ap_list:
        print '\nRadius ' + str(r_ap)

        t = Table()
        for i in xrange(len(ids)):
            photo_fn = os.path.expanduser(os.path.join(__root__,photom_dir,ids[i]+'_'+method+'_r'+str(r_ap)))
            
            photo = ascii.read(photo_fn, header_start=0, format='basic', delimiter=',' )
            
            if i==0:
                keys_info = ['Band', 'DATE', 'Frame', 'MJD', 'EXPOS', 'Airmass']
                name_info = ['#Band', 'Date', 'Frame', 'MJD', 'Exp_time', 'Airmass']
                for k in xrange(len(keys_info)):
                    t.add_column(Column(photo[keys_info[k]], name= name_info[k]))
        
            i_=str(i)
            keys_phot  = ['2MASS', 'X','Y', 'Flux', 'aperture_sum', 'aperture_badpix', 'aperture_area', 'background']
            name_phot = ['2massid','xpos','ypos','Flux','Rawflux_ap','Mask_pix_ap','Pix_ap','Bkgflux_ann']
            name_phot = [name+i_ for name in name_phot]
            for k in xrange(len(keys_phot)):
                t.add_column(Column(photo[keys_phot[k]], name= name_phot[k]))
        
        fs = os.path.expanduser(os.path.join(__root__,photom_dir, 
                                             'photometry_{:s}_rad{:s}_{:s}_v1_180115.'.format(band, str(round(r_ap,1)), method)+file_format))
        t.write(fs, format=file_format)