import nptipsyreader
import matplotlib.pyplot as plt
import numpy as np
import pdb

def rhot(tipsy):
    hist, ys, xs = np.histogram2d(np.log10(tipsy.gas['temp']), np.log10(tipsy.gas['rho']), bins=(100,100))
    xs = .5*(xs[:-1]+xs[1:])
    ys = .5*(ys[:-1]+ys[1:])
    scalemin = np.min(hist[hist>0])
    scalemax = np.max(hist[hist>0])
    levels = np.logspace(np.log10(scalemin), np.log10(scalemax), 50)
    plt.contour(xs, ys, hist, levels)
    plt.colorbar()
    plt.xlim(-3, 8)
    plt.ylim(2.5, 8)
    plt.xlabel('log(rho [cm-3])')
    plt.ylabel('log(temp [K])')
    plt.title(tipsy.filename + ' z ~ %0.1f' % (1./np.float(tipsy.t)-1))
    plt.savefig(tipsy.filename + '.rhot.png')
    plt.clf()

def metalstform(tipsy):
    recentsf = tipsy.star['tform'] >= (np.max(tipsy.star['tform']) - 0.05/tipsy.timeunit)
    sfr = np.sum(tipsy.star['mass'][recentsf])*np.float(tipsy.paramfile['dMsolUnit'])/5e7 #Msol/yr 
    hist,  ys, xs = np.histogram2d(tipsy.star['metals'], 
                                   tipsy.star['tform']*tipsy.timeunit, 
                                   bins=(100,100))
    xs = .5*(xs[:-1]+xs[1:])
    ys = .5*(ys[:-1]+ys[1:])
    scalemin = np.min(hist[hist>0])
    scalemax = np.max(hist[hist>0])
    levelbreak = (scalemax+scalemin)/2.
    levels = np.logspace(np.log10(scalemin), np.log10(scalemax), 25)
    #levels = np.concatenate((np.logspace(np.log10(scalemin), np.log10(levelbreak), 25), 
    #                         np.linspace(levelbreak, scalemax, 25)))
    plt.contour(xs, ys, hist, levels, label='csfr = %0.2f' % sfr)
    plt.colorbar()
    plt.legend()
    plt.xlabel('tform [Gyrs]')
    plt.ylabel('metals')
    plt.title(tipsy.filename + 'z ~ %0.1f' % (1./np.float(tipsy.t)-1))
    plt.savefig(tipsy.filename + '.metalstform.png')
    plt.clf()

def sfr(tipsy, tipsyindex):
    integrationtime = 5e7 #years
    recentt = np.max(tipsy.star['tform'])
    recentsf = tipsy.star['tform'][tipsyindex] >= (recentt - integrationtime/(1e9*tipsy.timeunit))
    sfr = np.sum(tipsy.star['mass'].astype('float64')[tipsyindex][recentsf])*np.float(tipsy.paramfile['dMsolUnit'])/integrationtime #Msol/yr
    #pdb.set_trace()
    return sfr, recentt*tipsy.timeunit, 1./tipsy.t-1.
