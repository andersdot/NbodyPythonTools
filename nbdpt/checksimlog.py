import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import pdb

class Logoutput(object):

    def __init__(self, *arg, **kwargs):
        self.filename = glob.glob('*.log')[0]
        self.title = ".".join(self.filename.split('.')[:-1])

        self.names = ('time', 'redshift', 'TotalEVir', 'TotalE', 'Kinetic', 'Virial', 'Potential', 'TotalECosmo', 'Ethermal', 'Lx', 'Ly', 'Lz', 'Wallclock')

    def read_param(self):
        self.paramfile = {}
        paramfilename = glob.glob('*.param')[0]
        f = open(paramfilename, 'rb')
        
        for line in f:
            try:
                if line[0] != '#' :
                    s = line.split('#')[0].split()
                    self.paramfile[s[0]] = " ".join(s[2:])
            except IndexError, ValueError:
                pass

        dKpcUnit = np.float(self.paramfile['dKpcUnit'])
        dMsolUnit = np.float(self.paramfile['dMsolUnit'])
        denunit = dMsolUnit/dKpcUnit**3
        self.vunit   = 8.0285*np.sqrt(6.6743e-8*denunit)*dKpcUnit
        self.h       = np.float(self.paramfile['dHubble0'])*10.*self.vunit/dKpcUnit
        self.timeunit=np.sqrt((dKpcUnit*3.086e21)**3/(6.67e-8*dMsolUnit*1.99e33))/(3600.*24.*365.24*1e9)



    def plot_log(self, plotfilename=None):
        time      = []
        redshift  = []
        fulltime = []
        fullredshift = []
        wallclock = []
        i = 0
        j = 0
        linenumber = []
        first = 'st'
        second = 'nd'
        third = 'rd'
        fourth = 'th'
        colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, 20)))
        
        for line in open(self.filename):
            entries = line.split()
            if (line[0] == '#'):
                if entries[1] == 'Running': 
                    ncpus = entries[3]
                    continue
                if (entries[1] == './ChaNGa.cosmo') or (entries[1] == './ChaNGa.Oct5'):
                    walltime_place = np.where(np.array(entries) == '-wall')[0]
                    walltime = entries[walltime_place + 1]
                if (time): #plot the stuff
                    fulltime.extend(time)
                    fullredshift.extend(redshift)
                    time      = np.array(time,      dtype='float')
                    redshift  = np.array(redshift,  dtype='float')
                    wallclock = np.array(wallclock, dtype='float')
                    linenumber = np.array(linenumber, dtype='float')
                    color = linenumber/np.max(linenumber)*255.
                    add_string = 'th'
                    if j == 1: add_string = first
                    if j == 2: add_string = second
                    if j == 3: add_string = third
                    print np.mean(wallclock)
                    if np.max(time*self.timeunit > 6.5): plt.scatter(time*self.timeunit, wallclock, color=next(colors), lw=0, label=np.str(j))
                    time = []
                    redshift = []
                    wallclock = []
                    linenumber = []
                    j = j + 1
                    continue
                if not time: continue
            time.append(entries[0])
            redshift.append(entries[1])
            wallclock.append(entries[12])
            linenumber.append(i)
            i = i+1

        fulltime.extend(time)
        fullredshift.extend(redshift)
        time      = np.array(time,      dtype='float')
        redshift  = np.array(redshift,  dtype='float')
        wallclock = np.array(wallclock, dtype='float')
        linenumber = np.array(linenumber, dtype='float')
        fulltime = np.array(fulltime, dtype='float')
        fullredshift = np.array(fullredshift, dtype='float')
        color = linenumber/np.max(linenumber)*255.


        plt.scatter(time*self.timeunit,wallclock, color = next(colors), s=80, lw=0, label=np.str(j))
        #plt.legend(loc=2)
        plt.xlabel('Time [Gyrs]')
        plt.ylabel('Wallclock Time')
        ax1 = plt.gca()
        ax2 = ax1.twiny()
        ind = np.linspace(0, len(fulltime)-1,10).astype(int)
        ax2.set_xticks(fulltime[ind]*self.timeunit)
        redlabel = ['%.1f' % fullredshift[i] for i in ind]
        ax2.set_xticklabels(redlabel)
        ax2.set_xlabel('redshift')
        plt.text(0.5, 1.085, self.title,
                 horizontalalignment='center',
             fontsize=20,
             transform = ax2.transAxes)
        ax1.set_xlim(0,fulltime[ind[-1]]*self.timeunit)
        if (plotfilename): plt.savefig(plotfilename+'.png')
        else: plt.show()
    
