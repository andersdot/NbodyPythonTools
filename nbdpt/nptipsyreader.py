import xdrlib
import numpy as np
import glob
import struct
import pdb

class Tipsy(object):
    def __init__(self, filename):
        self.filename = filename
        #self.f = open(filename)
        """ if self.paramfile.get('bDoublePos',0) : 
            ptype = 'd'
        else : 
            ptype = 'f'
            
        if self.paramfile.get('bDoubleVel',0) : 
            vtype = 'd'
        else : 
            vtype = 'f'"""
        ptype = '>f4'
        vtype = '>f4'
        ptype_swap = '<f4'
        vtype_swap = '<f4'
        self.byte_swap = False
        self._g_dtype = np.dtype({'names': ("mass","x","y","z","vx","vy","vz",
                                            "rho","temp","hsmooth","metals","phi"),
                                  'formats': ('>f4',ptype,ptype,ptype,vtype,vtype,
                                              vtype,'>f4','>f4','>f4','>f4','>f4')})
        self._d_dtype = np.dtype({'names': ("mass","x","y","z","vx","vy","vz",
                                            "eps","phi"),
                                  'formats': ('>f4',ptype,ptype,ptype,vtype,vtype,
                                              vtype,'>f4','>f4')})
        self._s_dtype = np.dtype({'names': ("mass","x","y","z","vx","vy","vz",
                                            "metals","tform","eps","phi"),
                                  'formats': ('>f4',ptype,ptype,ptype,vtype,vtype,
                                              vtype,'>f4','>f4','>f4','>f4')})

        self._g_dtype_swap = np.dtype({'names': ("mass","x","y","z","vx","vy","vz",
                                            "rho","temp","hsmooth","metals","phi"),
                                  'formats': ('<f4',ptype_swap,ptype_swap,ptype_swap,vtype_swap,vtype_swap,
                                              vtype_swap,'<f4','<f4','<f4','<f4','<f4')})
        self._d_dtype_swap = np.dtype({'names': ("mass","x","y","z","vx","vy","vz",
                                            "eps","phi"),
                                  'formats': ('<f4',ptype_swap,ptype_swap,ptype_swap,vtype_swap,vtype_swap,
                                              vtype_swap,'<f4','<f4')})
        self._s_dtype_swap = np.dtype({'names': ("mass","x","y","z","vx","vy","vz",
                                            "metals","tform","eps","phi"),
                                  'formats': ('<f4',ptype_swap,ptype_swap,ptype_swap,vtype_swap,vtype_swap,
                                              vtype_swap,'<f4','<f4','<f4','<f4')})

    def _read_param(self):
        try:
            paramfilename = [f for f in glob.glob('*.param') if re.match('^(cosmo|h)', f)]
        except IndexError:
            try:
                paramfilename = [f for f in glob.glob('../*.param')
 if re.match('^(cosmo|h)', f)_]               print 'There is no param file in this directory, trying one up'
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
        try: hub = np.float(paramfile['dHubble0'])
        except KeyError: hub=0.
        dunit = np.float(paramfile['dKpcUnit'])
        munit = np.float(paramfile['dMsolUnit'])
        denunit = munit/dunit**3.
        self.velunit = 8.0285*np.sqrt(6.6743e-8*denunit)*dunit
        hubunit = 10.*self.velunit/dunit
        self.h = hub*hubunit
        f.close()

    def unpack_header(self):
        f = open(self.filename)
        (t, nbodies, ndim, nsph, ndark, nstar) = struct.unpack('>diiiii', f.read(28))
        if (ndim < 1) or (ndim>3):
            f.seek(0)
            (t, nbodies, ndim, nsph, ndark, nstar) = struct.unpack('<diiiii', f.read(28))
            self.byte_swap=True
        f.close()
        return (t, nbodies, ndim, nsph, ndark, nstar)
        
    def unpack_gas(self, skip=None):
        if skip:
            ar = np.empty(self.ngas, dtype=self._g_dtype)
        else:
            f = open(self.filename)
            f.seek(struct.calcsize('>diiiiii'))
            if self.byte_swap:dtype = self._g_dtype_swap
            else:          dtype = self._g_dtype
            ar = np.core.records.fromfile(f, dtype=dtype, shape=self.ngas)
            f.close()
        return ar
        
    def unpack_dark(self, skip=None):
        if skip:
            ar = np.empty(self.ndark, dtype=self._d_dtype)
        else:
            f = open(self.filename)
            f.seek(struct.calcsize('>diiiiii')+struct.calcsize('12f4')*self.ngas)
            if self.byte_swap: dtype=self._d_dtype_swap
            else:           dtype=self._d_dtype
            ar = np.core.records.fromfile(f, dtype=dtype, shape=self.ndark)
            f.close()
        return ar

    def unpack_star(self, skip=None):
        if skip:
            ar = np.empty(self.nstar, dtype=self._s_dtype)
        else:
            f = open(self.filename)
            f.seek(struct.calcsize('>diiiiii')+struct.calcsize('12f4')*self.ngas+struct.calcsize('9f4')*self.ndark)
            if self.byte_swap: dtype=self._s_dtype_swap
            else:           dtype=self._s_dtype
            ar = np.core.records.fromfile(f, dtype=dtype, shape=self.nstar)
            f.close()
        return ar

    def _read(self, skipgas=None, skipdark=None, skipstar=None):
        (self.t, self.nbodies, self.ndim, self.ngas, self.ndark, self.nstar) = self.unpack_header()
        self.gas = self.unpack_gas(skip=skipgas)
        self.dark = self.unpack_dark(skip=skipdark)
        self.star = self.unpack_star(skip=skipstar)
