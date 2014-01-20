import pdb
#import struct
import numpy as np
#from astroML.plotting import scatter_contour
import matplotlib.pyplot as plt
#from scipy.integrate import quad
import cosmology
import pynbody
from .. import readparam
from .. import starlog as  sl
from numba.decorators import jit, autojit
import numba 
nftype = numba.double[:]
nitype = numba.int32[:]

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

def plot():
	volumesize = 6.25
	params = readparam.Params()
	cosmo = cosmology.Cosmology(float(params.paramfile['dOmega0']), params.h)
	cosmo.Age()
	starlog = sl.starlog()
	
	
	z = np.arange(0, 20, 0.05)
	times = np.empty_like(z)
	tform = starlog['tform']*params.timeunit
	for i in range(len(z)): times[i] = cosmo.AgeatZ(z[i])
	zformind = findz(times, tform)
	zform = z[zformind]
	dz = 1
        zbins = np.arange(np.min(zform), np.max(zform), dz)
        tbins = np.empty_like(zbins)
        for i in range(len(zbins)): tbins[i] = cosmo.AgeatZ(zbins[i])
        massform = starlog['massform']*float(params.paramfile['dMsolUnit'])
        boxsize = np.float(params.paramfile['dKpcUnit'])/1000.*(params.h/0.7)
        weights = massform/boxsize**3.
	hist, bins = np.histogram(zform, bins=zbins, 
				  weights=weights)
	x =  (bins[1:]+bins[:-1])/2.
        yerr = (x<0.9)*0.13 + ((x>=0.9)&(x<1.5))*0.17 + ((x>=1.5)&(x<3.))*0.19 + ((x>=3.)&(x<=8.))*.27 + (x>8.)*0.
        timenorm = (tbins[:-1] - tbins[1:])
        err = 1. + (hist + 0.75)**0.5
	plt.errorbar(x, np.log10(hist/(timenorm*1e9)), yerr=err/(timenorm*1e9), linewidth=2, drawstyle='steps-mid')
	zz = [4, 5, 6, 7]
	ss = [-1.14,-1.33,-1.58,-1.78]
	plt.scatter(zz, ss, color='g', label='Duncan, K et al 2014')
	A = -0.997
	B = 0.241
	C = 0.180
	z0 = 1.24
	plt.errorbar(x, (np.log10(csfr(x))), yerr=yerr, linestyle="None", marker="None", color='r') 
	plt.plot(z, (np.log10(csfr(z))), color='r', label='Behroozi et al 2013') 
	plt.legend()
	plt.xlabel('z')
	plt.ylabel('M$_{\odot}$ yr$^{-1}$ Mpc$^{-3}$')
        plt.fill([8,17,17,8], [-3.5,-3.5,-0.5,-0.5], 'grey', alpha=0.2, edgecolor='grey')
        plt.xlim(0, 17)
        plt.ylim(-3.5, -0.5)
        plt.show()
        pdb.set_trace()
if __name__ == '__main__': plot()
