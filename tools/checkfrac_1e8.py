import nptipsyreader
import numpy as np
import plots
import glob
import matplotlib.pyplot as plt
import pdb
import starlog
import cosmography
import fof2stat

threshold = 1e8

outputs = glob.glob('*.stat')
outputs.sort()

stats   = glob.glob('*.stat')
stats.sort()

fraction = np.empty(len(outputs), dtype='float')
ct = np.empty(len(outputs), dtype='float')

for i in range(len(outputs)):
	output = ".".join(outputs[i].split('.')[:-4])
	statdtype = np.dtype({'names': ('grp'  ,'mvir'   ,'npart'   ,'ngas'    ,'nstar',
					'ndark','gasmass','starmass','darkmass', 'rmax', 
					'x'    , 'y'     , 'z'      , 'vx'     , 'vy'  , 
					'vz'),
			      'formats': ('int32'  ,'float64','int32'  ,'int32'  ,'int32'  ,
					  'int32'  ,'float64','float64','float64','float64',
					  'float64','float64','float64','float64','float64',
					  'float64')})
	stat = np.genfromtxt(stats[i], dtype=statdtype, skip_header=1)
	tipsy = nptipsyreader.Tipsy(output)
	tipsy._read_param()
	tipsy._read(skipdark=True, skipgas=True)		
	ct[i] = np.max(tipsy.star['tform'])*tipsy.timeunit
	fraction[i] = np.sum(stat[stat['starmass'] > threshold]['nstar'])/np.float(tipsy.nstar)

plt.clf()
plt.plot(ct, fraction)
plt.xlabel('Time [Gyr]')
plt.ylabel('Fraction')
plt.show()
plt.savefig('fractions.png')
plt.clf()


