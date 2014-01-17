import pdb
#import struct
import numpy as np
#from astroML.plotting import scatter_contour
import matplotlib.pyplot as plt
#from scipy.integrate import quad
import cosmology
import pynbody
import nbdpt.readparam as readparam
import nbdpt.starlog as  starlog
from numba.decorators import jit, autojit

@jit(argtypes=(nftype, nftype))
def findz(times, tformstar):

     length = tformstar.shape[0]
     zformind = np.zeros(length, dtype=np.int32)

     for i in range(length): 
          zformind[i] = np.abs(times - tformstar[i]).argmin()

     return zformind


def csfr(z):
	z0 = 1.243
	A = -0.997
	B = 0.241
	C = 0.180
	return C/(10.**(A*(z-z0)) + 10.**(B*(z-z0)))

def plot(simfilename):
	volumesize = 6.25
	params = readparam.Params()
	sim = pynbody.load(filename)
	cosmo = cosmology.Cosmo(float(params.paramfile['dOmegaM']), params.h)
	cosmo.Age()
	starlog = starlog.starlog()
	
	
	z = np.arange(0, 20, 0.05)
	times = np.empty_like(z)
	tform = starlog['tform']*param.timeunit
	for i in range(len(z)): times[i] = cosmo.AgeatZ(z[i])
	zformind = findz(times, tform)
	zform = z[zformind]
	dz = 1
	hist, bins = np.histogram(zform, bins=(np.max(zform)-np.min(zform))/dz, 
				  weight = starlog['massform']*float(params.paramfile['dMsolUnit']))
	normalization = (bins[1:]-bins[:-1])*(np.float(params.paramfile['dKpcUnit'])/1000.
                                              *(params.h/0.7))**3
	x =  (bins[1:]+bins[:-1])/2.
        err = 1. + (hist + 0.75)**0.5
	plt.errorbar(x, np.log10(hist/normalization), linewidth=2, drawstyle='steps-mid')
	zz = [4, 5, 6, 7]
	ss = [-1.14,-1.33,-1.58,-1.78]
	plt.plot(zz, ss, 'Duncan, K et al 2014')
	A = -0.997
	B = 0.241
	C = 0.180
	z0 = 1.243
	plt.error
	plt.plot(z, (np.log10(csfr(z))), label='Behroozi') 
	plt.legend()
	plt.xlabel('z')
	plt.ylabel('M$_{\odot}$ yr$^{-1}$ Mpc$^{-3}$')
	plt.savefig('cosmo6BHC52_csfr.png')

if __name__ == '__main__': plot()
