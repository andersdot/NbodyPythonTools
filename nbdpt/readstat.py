import numpy as np

def readstat(filename, amiga=False):
    if amiga: statdtype = np.dtype({'names': ('grp'  ,'npart'  ,'ngas'    ,'nstar',
                                              'ndark','mvir','gasmass','starmass','darkmass', 'rmax',
                                              'x'    , 'y'     , 'z'      , 'vx'     , 'vy'  ,
                                              'vz'),
                                    'formats': ('int32'  ,'int32','int32'  ,'int32'  , 'int32'  ,
                                                'float64'  ,'float64','float64','float64', 'float64',
                                                'float64','float64','float64','float64', 'float64',
                                                'float64')})
    else: statdtype = np.dtype({'names':('grp', 'npart', 'ngas', 'nstar', 'ndark', 
                                         'mvir', 'rvir', 'gasmass', 'starmass', 'darkmass', 
                                         'vmax', 'rmax', 'sigmav', 'x', 'y', 'z', 'vx', 
                                         'vy', 'vz','contaminated', 'satellite', 'false'),
                                'formats':tuple(5*'int32 '.split() + 14*'float64 '.split() + ['|S10', 'bool', '|S10'])})
    return np.genfromtxt(filename, dtype=statdtype, skip_header=1)
