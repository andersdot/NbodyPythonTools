import numpy as np
import pdb
import matplotlib.pyplot as plt
import glob
from .. import nptipsyreader

class profile(object):
    def __init__(self, filename):
        self.filename = filename
        tipsyfile = glob.glob('*.den')[0]
        tipsy = nptipsyreader.Tipsy(tipsyfile)
        tipsy._read_param()
        self.length_unit = np.float(tipsy.paramfile['dKpcUnit'])
        self.mass_unit = np.float(tipsy.paramfile['dMsolUnit'])
    def read(self):
        f = open(self.filename)
        firstline = f.readline()
        print firstline
        print firstline.strip('\n')
        print firstline.strip('\n').split()
        self.names = firstline.strip('\n').split()[1:]
        f.close()
        
        profile = np.genfromtxt(self.filename, names=self.names, skiprows=1)
        return profile
    
if __name__=='__main__':

    starprofile = profile('cosmo6.25PLK.z.0.mainhalo.star.dat')
    starinfo    = starprofile.read()
    
    gasprofile = profile('cosmo6.25PLK.z.0.mainhalo.gas.dat')
    gasinfo    = gasprofile.read()
    
    dmprofile = profile('cosmo6.25PLK.z.0.mainhalo.dark.dat')
    dminfo    = dmprofile.read()
    
    allprofile = profile('cosmo6.25PLK.z.0.mainhalo.all.dat')
    allinfo   = allprofile.read()
    
    x = 'radius'
    y = 'rho'

    plt.plot(allinfo[x]*allprofile.length_unit, allinfo[y]*allprofile.mass_unit/allprofile.length_unit**3., label='all')
    plt.plot(dminfo[x]*dmprofile.length_unit, dminfo[y]*dmprofile.mass_unit/dmprofile.length_unit**3., label='dark')
    plt.plot(gasinfo[x]*gasprofile.length_unit, gasinfo[y]*gasprofile.mass_unit/gasprofile.length_unit**3., label='gas')
    plt.plot(starinfo[x]*starprofile.length_unit, starinfo[y]*starprofile.mass_unit/starprofile.length_unit**3., label='star')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('radius [kpc]')
    plt.ylabel('$rho$ [M$_{\odot}$ kpc$^{-3}$]')
    plt.title('Cosmo6.25PLK Main Halo Density Profile')
    plt.xlim(1e-2, 20)
    plt.show()
    pdb.set_trace()
    
