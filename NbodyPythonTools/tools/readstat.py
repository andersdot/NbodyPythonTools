import numpy as np

def readstat(filename):
    statdtype = np.dtype({'names': ('grp'  ,'mvir'   ,'npart'   ,'ngas'    ,'nstar',
                                    'ndark','gasmass','starmass','darkmass', 'rmax',
                                    'x'    , 'y'     , 'z'      , 'vx'     , 'vy'  ,
                                    'vz'),
                          'formats': ('int32'  ,'float64','int32'  ,'int32'  , 'int32'  ,
                                      'int32'  ,'float64','float64','float64', 'float64',
                                      'float64','float64','float64','float64', 'float64',
                                      'float64')})
    
    
    return np.genfromtxt(filename, dtype=statdtype, skip_header=1)
