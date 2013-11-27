import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import nptipsyreader
import numpy as np
import pdb

statfiles = glob.glob('*.stat')
statfiles.sort()

statdtype = np.dtype({'names': ('grp'  ,'mvir'   ,'npart'   ,'ngas'    ,'nstar',
                                'ndark','gasmass','starmass','darkmass', 'rmax',
                                'x'    , 'y'     , 'z'      , 'vx'     , 'vy'  ,
                                'vz'),
                      'formats': ('int32'  ,'float64','int32'  ,'int32'  , 'int32'  ,
                                  'int32'  ,'float64','float64','float64', 'float64',
                                  'float64','float64','float64','float64', 'float64',
                                  'float64')})

datadir = '/u/sciteam/lmanders/data/release-sfh_z0_z8_052913/smmr/'
dataprefix = 'c_smmr_z'
datasuffix = '0_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'

colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, 11))) #len(statfiles))))

for i in range(len(statfiles)):
    stat = np.genfromtxt(statfiles[i], dtype=statdtype, skip_header=1)
    tipsyfile = ('.').join(statfiles[i].split('.')[0:4])
    tipsy = nptipsyreader.Tipsy(tipsyfile)
    tipsy._read(skipdark=True, skipgas=True)
    tipsy._read_param()
    time = np.max(tipsy.star['tform'])*tipsy.timeunit
    redshift = 1./tipsy.t -1.
    if np.round(redshift, 1) in [0.1, 1, 2, 3, 4, 5, 6, 7, 8]:

        datafile = datadir+dataprefix + str(np.round(redshift)) + datasuffix
        (mvir, smhm, blah, blah) = np.genfromtxt(datafile, unpack=True)
        nonzero = (stat['starmass'] > 0.) & (stat['starmass'] > 4.*np.float(tipsy.paramfile['dInitStarMass'])*np.float(tipsy.paramfile['dMsolUnit']))
        thiscolor = next(colors)
        mass = stat['mvir'][nonzero]*(0.6777/0.7)*1.2
        starm = stat['starmass'][nonzero]*(0.6777/0.7)*0.7
        plt.plot(10.**mvir, 10.**smhm*10.**mvir, color = thiscolor)
        plt.scatter(mass, starm, color=thiscolor, lw=0, label='z~' + "{:.2}".format(redshift)) #"{:.4}".format(time) + '[Gyrs] z~' + "{:.2}".format(redshift))

plt.yscale('log')
plt.xscale('log')
plt.xlabel('M$_{h}$ [M$_{\odot}$]', fontsize=15)
plt.ylabel('M$_{*}$/M$_{h}$', fontsize=15)
plt.ylim(1e5,)
plt.title('Cosmo6.25PLK Evolution of Stellar Mass Fractions')
plt.legend()
plt.show()
