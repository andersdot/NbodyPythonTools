import pynbody
import numpy as np
import matplotlib.pyplot as plt
import glob

if __name__ == '__main__':
    filenames = glob.glob('*.grp')
    
    for i in range(len(filenames)): 
        filename = '.'.join(filenames[i].split('.')[0:3])
        sim = pynbody.load(filename)
        halos = sim.halos()
        sim.stars['massform'].units = sim._paramfile['dMsolUnit'] + ' Msol'
        sim.physical_units()
        
        nhalos = np.max(sim['amiga.grp'])
        halo_mag = np.zeros(nhalos)
        halo_mass = np.zeros(nhalos)
        halo_lum = np.zeros(nhalos)
        
        for i in range(nhalos-1): 
            if len(halos[i+1].star) > 0:
                halo_mag[i] = pynbody.analysis.luminosity.halo_mag(halos[i+1], band='u')
                halo_mass[i] = np.sum(halos[i]['mass'])
                halo_lum[i] = pynbody.analysis.luminosity.halo_lum(halos[i+1], band='u')
                
        np.savetxt(filename + '.mag', np.column_stack((halo_mass, halo_lum, halo_mag)), header= 'Mass [Msol]   Luminosity [Lsol]    AbsoluteMagnitude', fmt='%.5e')


"""
plt.plot(halo_mass, halo_mag, 'ko')
plt.show()
plt.plot(halo_mass, halo_lum, 'bo')
plt.show()
plt.plot(halo_lu, halo_mag, 'go')
plt.show()

comovingboxsize = sim.properties['boxsize'].in_units('Mpc')/sim.properties['a']

#z ~ 2.9, lambda ~ 1700
hist, bins = np.histogram(halo_mag[~np.isnan(halo_mag)],bins=30)
normalization = (bins[1:]-bins[:-1])*comovingboxsize**3
err = 1. + (hist + 0.75)**0.5
plt.errorbar((bins[1:]+bins[:-1])/2., hist/normalization, yerr = err/normalization, label='Cosmo50 z~3.5')
phistar = 1.06e-3
mstar = -21.11
alpha = -1.94

phi = 2./5*phistar*np.log(10)*(10.**(2./5.*(mstar - bins)))**(alpha+1.)*np.exp(-10.**(2./5.*(mstar-bins)))

x = np.array((-23.5, -22.9, -22.3, -21.7, -20.5, -19.99, -19.4, -18.8, -18.25, -17.7, -17.1, -16.5))
y = np.array((-5.1, -4.05, -3.7, -3.29, -2.85, -2.61, -2.35, -2.21, -2.0, -2.0, -1.69, -1.6))

plt.plot(bins, phi, label='z~3 Bian 2013')
plt.yscale('log')
plt.xlim(-10, -26)
plt.ylim(1e-5, 1)
plt.plot(x, 10.**y, label = 'z~4 Dust Corrected, Smit 2012')
plt.legend()
plt.xlabel('M$_{UV,AB}$')
plt.ylabel('log$_{10}$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.show()



np.savetxt(filename + '.mag', (halo_mass, halo_lum, halo_mag), header= 'Mass [Msol]   Luminosity [Lsol]    AbsoluteMagnitude', fmt='%.5e')
"""
