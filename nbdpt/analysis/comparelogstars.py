### Friday fix ahf->grp/stat
### Make more Salsa outputs
### Finish this program

import numpy as np
import os
import glob
import struct
import matplotlib.pyplot as plt

if __name__ == '__main__':
#### Read in param file to get units ####
    x = os.path.curdir
    parafile = glob.glob(os.path.join(x, '*.param'))[0]
    f = open(filename)
### convert simtime to real time ###
### add redshift axis ##############
    
##### Read in starlog file #########
    filename = glob.glob(os.path.join(x, '*.starlog'))[0]
    f = open(filename)
    
    file_structure = np.dtype({'names': ("iord","iorderGas","tform",
                                         "x","y","z",
                                         "vx","vy","vz",
                                         "massform","rhoform","tempform"),
                               'formats':('i4','i4','f8',
                                          'f8','f8','f8',
                                          'f8','f8','f8',
                                          'f8','f8','f8')})
    
    size = struct.unpack(">i", f.read(4))
    datasize = os.path.getsize(filename)-f.tell()
    g = np.fromstring(f.read(datasize), dtype=file_structure).byteswap()
    
    
    tmp = g['iord'].flatten()
    perm = tmp.argsort()
    aux = tmp[perm]
    flag = np.concatenate(([True],aux[1:]!=aux[:-1]))
    iord = aux[flag]
    indices = perm[flag]
    
    timeunit=np.sqrt((dKpcUnit*3.086e21)**3/(6.67e-8*dMsolUnit*1.99e33))/(3600.*24.*365.24)
    binnorm = 1. # 1./width of bin in years 1/actual normalization since put in weights
    tforms = g['tform'][indices]
    weight = g['massform'][indices]*binnorm
    sfhist, thebins, patches = plt.hist(tforms, weights=weigh, bins=bins, histtype='step')
    
    
#####################################
# plot the two on top of each other#
    
