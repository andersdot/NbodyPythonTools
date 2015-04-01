import glob
import struct
import os
import numpy as np
import pdb

def starlog():
    starlogfile = glob.glob('*.starlog')[0]
    f = open(starlogfile)
    
    size = struct.unpack('>i', f.read(4))[0]
    
    file_structure =  np.dtype({'names': ("iord","iorderGas","tform",
                                          "x","y","z",
                                          "vx","vy","vz",
                                          "massform","rhoform","tempform"),
                                'formats':('>i4','>i4','>f8',
                                           '>f8','>f8','>f8',
                                           '>f8','>f8','>f8',
                                           '>f8','>f8','>f8')})
    
    if size != file_structure.itemsize:
        print 'Starlog has wrong format'
        return
    
    shape = os.path.getsize(starlogfile)/size
    ar = np.core.records.fromfile(f, dtype=file_structure, shape=shape)
    """
    tmp = ar['iord'].flatten()
    perm = tmp.argsort()
    aux = tmp[perm]
    flag = np.concatenate(([True],aux[1:]!=aux[:-1]))
    iord = aux[flag]
    indices = perm[flag]
    ar = ar[indices]
    

    print len(ar)
    print len(np.unique(ar['iord']))
    """
    return ar
