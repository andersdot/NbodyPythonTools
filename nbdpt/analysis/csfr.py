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
from .. import findtipsy
from .. import readstat
from .. import nptipsyreader

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

def correctedsfr(threshold=1e8):
     tipsyfile = findtipsy.find()
     tipsyfile.sort()
     currentsfrfix = np.empty(len(tipsyfile), dtype='float64')
     currentsfr = np.empty_like(currentsfrfix)
     ct = np.empty_like(currentsfrfix)
     z = np.empty_like(currentsfrfix)
     fraction = np.empty_like(currentsfrfix)
     fractionstars = np.empty_like(currentsfrfix)

     for i in range(len(tipsyfile)):
          fgrp = open(tipsyfile[i] + '.rockstar.grp')
          print fgrp.readline()
          grp = [line.split('\n')[0] for line in fgrp]
          fgrp.close()
          grp = np.array(grp, dtype='int32')
          print grp[0]
          stat = readstat.readstat(tipsyfile[i] + '.rockstar.stat')
          wantedgrp = stat['grp'][stat['starmass']>threshold]
          
          tipsy = nptipsyreader.Tipsy(tipsyfile[i])
          tipsy._read()
          tipsy._read_param()
          z[i] = 1./tipsy.t -1.
          starindex = []
          for j in wantedgrp:
               tipsyindex = np.where(grp == j)[0]
               fullindex = np.where(tipsyindex > (tipsy.ngas + tipsy.ndark-1))[0]
               starindex.extend(tipsyindex[fullindex] - (tipsy.ngas + tipsy.ndark))
          starindex = np.array(starindex, dtype='int32')
          integrationtime = 5e7
          ct[i] = np.max(tipsy.star['tform'])*tipsy.timeunit
          recentsf = tipsy.star['tform'][starindex]*tipsy.timeunit >= (ct[i] - integrationtime/1e9)
          currentsfrfix[i] = np.sum(tipsy.star['mass'].astype('float64')[starindex][recentsf])*np.float(tipsy.paramfile['dMsolUnit'])/integrationtime
          
          recentsf = tipsy.star['tform']*tipsy.timeunit >= (ct[i] - integrationtime/1e9)
          currentsfr[i] = (np.sum(tipsy.star['mass'].astype('float64')[recentsf])*
                           np.float(tipsy.paramfile['dMsolUnit'])/integrationtime)
          fractionstars[i] = np.sum(stat[stat['starmass']>threshold]['nstar'])/float(np.sum(stat['nstar']))
          pdb.set_trace()
     return currentsfrfix, currentsfr, ct, z, fractionstars
     

def plot():

     sfrfix, sfr, ct, z, fstar = correctedsfr()
     plt.plot(z, sfrfix, label='fixed')
     plt.plot(z, sfr, label='unfixed')
     plt.legend()
     plt.yscale('log')
     plt.show()
     plt.plot(z, fstar, label='fractionstar')
     plt.legend()
     plt.show()
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
     plt.errorbar(x, np.log10(hist/(timenorm*1e9)), yerr=err/(timenorm*1e9), color='k', linewidth=2, drawstyle='steps-mid')
     zz = [4, 5, 6, 7]
     ss = [-1.14,-1.33,-1.58,-1.78]
     plt.errorbar(zz, ss, yerr=np.array((0.32, 0.16, 0.2, 0.14)), color='g', linestyle='None', marker='o', label='Duncan, K et al 2014')
     plt.errorbar(2., np.log10(0.148), yerr=[[0.063], [0.0627]], color='b', linestyle='None', marker='o',label='Alavi et al 2012')
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
     
