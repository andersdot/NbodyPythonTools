from .. import nptipsyreader
import numpy as np
import plots
import glob
import matplotlib.pyplot as plt
import pdb
from .. import starlog
from .. import cosmography
from ..fof import fof2stat

def csfr(z):
	z0 = 1.243
	A = -0.997
	B = 0.241
	C = 0.180
	return C/(10.**(A*(z-z0)) + 10.**(B*(z-z0)))

def outputsfr(threshold=1e8):
	outputs = glob.glob('*.grp')
	outputs.sort()
	grps    = glob.glob('*.grp')
	grps.sort()
	stats   = glob.glob('*.stat')
	stats.sort()
	currentsfrfix = np.empty(len(outputs), dtype='float')
	currentsfr = np.empty(len(outputs), dtype='float')
	ct   = np.empty(len(outputs), dtype='float')
	z    = np.empty(len(outputs), dtype='float')
	fraction = np.empty(len(outputs), dtype='float')
	fractionstars = np.empty(len(outputs), dtype='float')
	for i in range(len(outputs)):
		#read in grp file
		output = ".".join(outputs[i].split('.')[:-4])
		print output, grps[i], stats[i]
		fgrp = open(grps[i])
		grp = []
		grp = [line.split('\n')[0] for line in fgrp]
		fgrp.close()
		grp = np.array(grp, dtype='int32')
		grp = grp[1:]
		
		#read in stat file
		statdtype = np.dtype({'names': ('grp'  ,'mvir'   ,'npart'   ,'ngas'    ,'nstar',
						'ndark','gasmass','starmass','darkmass', 'rmax', 
						'x'    , 'y'     , 'z'      , 'vx'     , 'vy'  , 
						'vz'),
				      'formats': ('int32'  ,'float64','int32'  ,'int32'  ,'int32'  ,
						  'int32'  ,'float64','float64','float64','float64',
						  'float64','float64','float64','float64','float64',
						  'float64')})
		stat = np.genfromtxt(stats[i], dtype=statdtype, skip_header=1)
		wantedgrp = stat['grp'][stat['starmass'] > threshold]
		#read in tipsy file
		tipsy = nptipsyreader.Tipsy(output)
		tipsy._read_param()
		tipsy._read()		
		starindex = []

		for j in wantedgrp: 
			tipsyindex = np.where(grp == j)[0]
			fullindex = np.where(tipsyindex > (tipsy.ngas + tipsy.ndark- 1))[0]
			starindex.extend(tipsyindex[fullindex] - (tipsy.ngas + tipsy.ndark))
		starindex = np.array(starindex, dtype='int32')
		currentsfrfix[i], ct[i], z[i] = plots.sfr(tipsy,starindex)
		recentsf = tipsy.star['tform'] >= (np.max(tipsy.star['tform']) - 5e7/(1e9*tipsy.timeunit))
		currentsfr[i] = (np.sum(tipsy.star['mass'].astype('float64')[recentsf])*np.float(tipsy.paramfile['dMsolUnit'])/5e7)
		fractionstars[i] = np.sum(stat[stat['starmass']>1e8]['nstar'])/np.sum(stat['nstar'])
	
	sortarg = np.argsort(ct)
	ct = ct[sortarg]
	currentsfr = currentsfr[sortarg]
	z = z[sortarg]
	plt.clf()
	plt.plot(ct, fraction[sortarg], label='Fraction SFR')
	plt.plot(ct, fractionstars[sortarg], label='Fraction Stars')
	plt.xlabel('Time [Gyr]')
	plt.ylabel('Fraction')
	plt.savefig('fractions.png')
	plt.clf()
	return currentsfr, currentsfrfix, ct, z, tipsy

if __name__ == '__main__':

### SFR from outputs ####
	currentsfr, currentsfrfix, ct, z, tipsy = outputsfr()
	normalization = (np.float(tipsy.paramfile['dKpcUnit'])/1e3*(tipsy.h/0.7))**3.
	plt.plot(ct, currentsfr/normalization, 'ko', label='Outputs')
	plt.plot(ct, currentsfrfix/normalization, 'bo', label='Fixed Outputs')
########################
	
### SFR from starlog ###
	starlog = starlog.starlog()
	formationtime = starlog['tform']*tipsy.timeunit #in Gyrs
	formationmass = starlog['massform']*np.float(tipsy.paramfile['dMsolUnit']) #in Msol
	binsize = 0.05 #50 Million years
	bins = np.arange(np.min(formationtime), np.max(formationtime), binsize)
	boxlength = np.float(tipsy.paramfile['dKpcUnit'])/1e3
	h         = tipsy.h
	binnorm = 1./(binsize*1e9*(boxlength*(h/0.7))**3.) #in years
	weight = formationmass*binnorm
	
	sfhist, thebins, patches = plt.hist(formationtime, weights=weight, bins=bins, 
					    histtype='step', label='Starlog')
#########################
	
### SFR from Behroozi ###
	redshift = np.arange(2000)/100.
	OmegaM = np.float(tipsy.paramfile['dOmega0'])
	OmegaL = np.float(tipsy.paramfile['dLambda'])
	cosmo = cosmography.Cosmology(OmegaM, OmegaL, OmegaM+OmegaL, tipsy.h)
	tt = 13.7 - cosmo.Lookback_Time(redshift)*13.7
	plt.plot(tt, csfr(redshift), label='Behroozi')
#plt.plot(np.zeros(1e6)+ 0.66, np.log10(np.arange(1e6)/1e6))
#########################
	
	plt.ylim(0, 1.2*np.max(sfhist))
	plt.xlim(0, np.max(formationtime))
	plt.xlabel('Time [Gyrs]', fontsize='large')
	plt.ylabel('SFR [M$_\odot$ yr$^{-1}$ Mpc$^{-3}$]', fontsize='large')
	plt.legend()
	
	
	x0, x1 = plt.gca().get_xlim()
	pz = plt.twiny()
	ind = np.arange(0, len(z), 2)
	pz.set_xticks(ct[ind]) #[ind])
	pz.set_xticklabels(['%.1f' % a for a in z[ind]])
	pz.set_xlabel('z')
	pz.set_xlim(x0, x1)
	plt.text(0.5, 1.085, 'Cosmic SFR',
		 horizontalalignment='center',
		 fontsize=20,
		 transform = pz.transAxes)
	
	
	
	plt.show()
