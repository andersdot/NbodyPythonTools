import xdrlib
import numpy as np
import glob
import struct
import pdb

class Params(object):    
    def __init__(self, *arg, **kwargs):
        try:
            paramfilename = glob.glob('*.param')[0]
        except IndexError:
            try:
                paramfilename = glob.glob('../*.param')[0]
                print 'There is no param file in this directory, trying one up'
            except IndexError:
                print "Can't find param file"
                return
        f = open(paramfilename, 'rb')
        paramfile = {}
        for line in f:
            try:
                if line[0] != '#':
                    s = line.split('#')[0].split()
                    paramfile[s[0]] = "".join(s[2:])
            except IndexError, ValueError:
                pass
        self.paramfile = paramfile
        dKpcUnit = np.float(paramfile['dKpcUnit'])
        dMsolUnit = np.float(paramfile['dMsolUnit'])
        self.timeunit=np.sqrt((dKpcUnit*3.086e21)**3/
                              (6.67e-8*dMsolUnit*1.99e33)
                              )/(3600.*24.*365.24*1e9)
        hub = np.float(paramfile['dHubble0'])
        dunit = np.float(paramfile['dKpcUnit'])
        munit = np.float(paramfile['dMsolUnit'])
        denunit = munit/dunit**3.
        self.velunit = 8.0285*np.sqrt(6.6743e-8*denunit)*dunit
        hubunit = 10.*self.velunit/dunit
        self.h = hub*hubunit
        f.close()
        
