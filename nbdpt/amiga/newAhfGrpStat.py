import numpy as np
import pynbody
import glob
import pdb
import struct
from .. import nptipsyreader

class AHF(object):
    def __init__(self, tipsyfile):
        self._filename = tipsyfile

    def _readtipsy(self):
        self.tipsy = tipsyreader.Tipsy(self._filename)
        self.tipsy._read_param()
        self.tipsy._read()

        self.boxsize = np.float(self.tipsy.paramfile['dKpcUnit'])
        self.munit   = np.float(self.tipsy.paramfile['dMsolUnit'])
        self.denunit = self.munit/self.boxsize**3
        self.vunit   = 8.0285*np.sqrt(6.6743e-8*self.denunit)*self.boxsize
        self.timeunit = self.boxsize/self.vunit*0.7781311
        self.h       = np.float(self.tipsy.paramfile['dHubble0'])*10.*self.vunit/self.boxsize
     
    def halos(self):
        # read halofile into memory
        tipsyfile = self._filename
        halofile = glob.glob(tipsyfile + '*halos')[0]
        f = open(halofile, 'r')
        haloColumnNames = f.readline().strip('#\n').split('\t')[:-1]
        int = 'i8'
        flt = 'f8'
        str = 'a8'
        
        dtype = (int,int,int,flt,int,flt,flt,flt,flt,flt,flt,flt,
                 flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,
                 flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,
                 flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,
                 flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,
                 flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,
                 flt,flt,flt,flt,flt,flt,flt,flt,flt,flt,flt)
        usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,16,18,19,20,35)
        
        halos = np.genfromtxt(halofile, names=haloColumnNames, dtype=dtype,usecols=usecols, skiprows=1)
        return halos

    def particles(self):
        #read particle file into memory
        tipsyfile = self._filename
        particlefile = glob.glob(tipsyfile + '*particles')[0]
        index = []
        for line in open(particlefile) : index.append(line.split()[0])
        index = np.array(index, dtype=np.long)
        return index

    def reorder(self, halos, index):
        #Put particles in gas, dark, star order
        particles = np.empty_like(index)
        grp       = np.empty_like(index)
        if self.tipsy.ngas != 0:
            darkindex = particles < self.tipsy.ndark
            nd        = np.sum(darkindex)
            
            starindex = (particles >= self.tipsy.ndark) & (particles < self.tipsy.ndark+self.tipsy.nstar)
            ns        = np.sum(starindex)
            
            gasindex  = particles >= self.tipsy.ndark+self.tipsy.nstar
            ng        = np.sum(gasindex)
            
            sparticles = np.empty(nd+ns+ng, dtype=np.long)
            sgrp       = np.empty(nd+ns+ng, dtype=np.long)
            
            sparticles[0:ng-1] = particles[gasindex] - self.tipsy.ndark - self.tipsy.nstar
            sgrp      [0:ng-1] = grp      [gasindex]
            
            sparticles[ng:ng+nd-1] = particles[darkindex] + self.tipsy.ngas
            sgrp      [ng:ng+nd-1] = grp      [darkindex]
            
            if nstar != 0:
                sparticles[nd+ng:nd+ng+ns-1] = particles[starindex]+ngas
                sgrp      [nd+ng:nd+ng+ns-1] = grp[starindex]
        else:
            sparticles = particles
            sgrp       = grp
        return sparticles, sgrp

    def grp(self, sparticles, sgrp):
        #Put the particle indices in ascending order, and get the grp number that each particle belongs to.
        sortind  = np.argsort(sparticles)
        sortgrp  = sgrp[sortind]
        totarray = np.empty(self.tipsy.nbodies, dtype=np.long)
        ind      = sparticles[sortind]
        for i in range(len(sortind)): totarray[ind[i]] = sortgrp[i]
        return totarray #grp file

    def write_grp(self, grp):
        #Write grp file to disk
        """NOTE: some particles are assigned by AMIGA to more than one halo!  This
        means that some of the ind values that go into totarray get written over
        after an intial assignment of grp number.  This is ok if the halos are
        read in by order of mass, as that guarantees that the reassigned value
        is a substructure grp number.  However, AMIGA currently orders the halos
        by npart rather than my mass.  However, it is likely that the oddball
        halos that are out of mass order are lowres halos that accidentally made
        it into the catalog.  So, for now I will not rearrange the input by mass.
        This may become necessary at a later time, though."""

        filename = self._filename+'.amiga.grp'
        f = open(filename, 'w')
        f.write(self.tipsy.nbodies)
        for i in range(len(grp)): f.write(grp[i])
        f.close()

    def stat(self, halos, sparticles, sgrp):
        sortind  = np.argsort(sparticles)
        sortgrp  = sgrp[sortind]
        sort     = sortgrp[np.argsort(sortgrp)]
        unique   = sort[np.uniq(sort)]
        ngrps    = len(unique)
        amigaind = unique - 1

        dtype = np.dtype({'names': ('grp','sat', 'contam', 'false', 'npart', 'ngas', 
                                    'nstar', 'ndar', 'mvir', 'rvir', 'gasmass', 'starmass', 
                                    'darkmass', 'vmax', 'rmax', 'sigv', 'xc', 'yc', 'zc', 'vx', 
                                    'vy', 'vz', 'aimgaorigid'),
                          'formats': ('l','a','a','a','l','l','l','l','d','d', 
                                      'd','d','d','d','d','d','d','d','d','d' 
                                      'd','l')})
        stat          = np.empty(ngrps, dtype = dtype) 
        stat['grp']   = unique
        stat['npart'] = halos['npart5'].npart[amigaind]
        stat['mvir']  = halos['mvir'][amigaind]/self.h
        stat['rvir']  = halos['rvir'][amigaind]/self.h
        stat['vmax']  = halos['vmax'][amigaind]
        stat['rmax']  = halos['rmax'][amigaind]/self.h
        stat['sigv']  = halos['sigv'][amigaind]
        stat['x']     = halos['x']   [aimgaind]/self.h
        stat['y']     = halos['y']   [aimgaind]/self.h
        stat['z']     = halos['z']   [aimgaind]/self.h
        stat['vx']    = halos['vx']  [amigaind]
        stat['vy']    = halos['vy']  [aimgaind]
        stat['vz']    = halos['vz']  [aimgaind]

        if self.tipsy.ngas+self.tipsy.nstar > 0:
            gmass = tipsy.gas['mass'][0]
            

if __name__ == 'main':
    tipsyfile = 'cosmo50p.512.016384'
    ahf   = AHF(tipsyfile)
    units = ahf.simunits()
    halos = ahf.halos()
    index = ahf.particles()
    
    boxsize = ahf.boxsize
    munit   = ahf.munit
    h0      = ahf.h
    vunit   = ahf.vunit
    ndark   = ahf.ndark
    ngas    = ahf.ngas
    nstars  = ahf.nstar

    particles = np.empty_like(index)
    grp       = np.empty_like(index)

    sortedparticles, sortedgrp = ahf.reorder(halos, index)
    grp = ahf.grp(sortedparticles, sortedgrp)
    ahf.write_grp(grp)


    
#define the initial and final index for each halo in the particle file
#the particle file is formatted such that the very first value is the total 
#number of halos in the file. For each halo, the number of partilces in the 
#halo is listed, and that many next lines contain the index of th eparticles 
#that belong to that halo. 
    """
    npart    = halos['npart5']
    startind = npart.cumsum() + np.arange(len(npart)) + 2
    np.insert(startind,0,1) #need to prepend the 0th index to the beginning of startind
    endind   = npart.cumsum() + np.arange(len(npart))
    
    mindm  = np.min(sim.dm['mass'])
    
    # get rid of the massive dm particles
    for i in range(len(startind)):
    halo = index[startind[i]:endind[i]]
    dark_ind = halo < ndark
    if dark_ind:
    dark_mass = np.sum(sim.dm[halo[dark_ind]]['mass'])
    massives  = sim.dm[halo[dark_ind]]['mass'] > 1.5*mindm
    if massives:
    heavy_dark_mass = np.sum(sim.dm[halo[dark_ind[massives]]]['mass'])
    if have_dark_mass/dark_mass > 0.1: continue
    else: particles[startind[i]:endind[i]] = index[startin[i]:endind[i]]
    else: particles[startind[i]:endind[i]] = index[startin[i]:endind[i]]
    else: continue
    grp[startind[i]:endind[i]] = i+1
    nonzero   = (particles != 0)
    particles = particles[nonzero]
    grp       = grp[nonzero]
    
    """
    
        
# Put particles indices in ascending order, and get grp number each particle belongs to

    pdb.set_trace()

"""

; **************************************************************************
; Now, create the .stat file
; **************************************************************************
sorted = sortedgrp(sort(sortedgrp))
unique = sorted(uniq(sorted))
ngrps = n_elements(unique)
amigaind = unique-1

stat = replicate({grp:0L, sat:strarr(1), contam:strarr(1), false:strarr(1), n_part:0L, n_gas:0L, n_star:0L, n_dark:0L, m_vir:dblarr(1), r_vir:dblarr(1), gas_mass:dblarr(1), star_mass:dblarr(1), dark_mass:dblarr(1), v_max:dblarr(1), r_max:dblarr(1), sig_v:dblarr(1), x_c:dblarr(1), y_c:dblarr(1), z_c:dblarr(1), vx_c:dblarr(1), vy_c:dblarr(1), vz_c:dblarr(1), amiga_orig_id:0L}, ngrps)

; Take care of the easy stuff first. 
; By definition, there should not be grp 0 
;for i=0L,ngrps-1 do stat[i].grp = i 
;for i=0L,ngrps-1 do stat[i].amiga_orig_id = amigaind[i]
for i=0L,ngrps-1 do stat[i].grp = unique[i] 
stat.n_part = npart(amigaind)
for i=0L,ngrps-1 do stat[i].m_vir = double(mvir(amigaind[i])/h0) 
for i=0L,ngrps-1 do stat[i].r_vir = double(rvir(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].v_max = double(vmax(amigaind[i]))
for i=0L,ngrps-1 do stat[i].r_max = double(rmax(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].sig_v = double(sigv(amigaind[i]))
for i=0L,ngrps-1 do stat[i].x_c = double(xc(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].y_c = double(yc(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].z_c = double(zc(amigaind[i])/h0)
for i=0L,ngrps-1 do stat[i].vx_c = double(vxc(amigaind[i]))
for i=0L,ngrps-1 do stat[i].vy_c = double(vyc(amigaind[i]))
for i=0L,ngrps-1 do stat[i].vz_c = double(vzc(amigaind[i]))
"""


"""
; Now figure out info on gas and stars, if they exist!
IF h.ngas+h.nstar GT 0 THEN BEGIN 
  gmass = double(g.mass[0])
  if h.nstar GT 0 then smass = double(s.mass[0])
  dmass = double(d.mass[0])
  FOR i=0L, ngrps-1 do begin
    particles_considered = index[startind(amigaind[i]):endind(amigaind[i])]
    gas_ind = where(particles_considered GT h.ndark+h.nstar)
    if gas_ind[0] NE -1 then stat[i].n_gas = n_elements(gas_ind) else stat[i].n_gas = 0
    if gas_ind[0] NE -1 then stat[i].gas_mass = total(gmass[particles_considered(gas_ind)-h.ndark-h.nstar]) else stat[i].gas_mass = 0
  IF h.nstar GT 0 then begin
    star_ind = where(particles_considered GE h.ndark and particles_considered LT h.ndark+h.nstar)
    if star_ind[0] NE -1 then stat[i].n_star = n_elements(star_ind) else stat[i].n_star = 0 
    if star_ind[0] NE -1 then stat[i].star_mass = total(smass[particles_considered(star_ind)-h.ndark]) else stat[i].star_mass = 0 
  ENDIF else stat[i].star_mass = 0
    dark_ind = where(particles_considered LT h.ndark)
    if dark_ind[0] NE -1 then begin
	stat[i].n_dark = n_elements(dark_ind)
	stat[i].dark_mass = total(dmass[particles_considered(dark_ind)])
	test2 = where(dmass[particles_considered(dark_ind)] GT mindm)
	if test2[0] NE -1 then stat[i].contam = 'contam' else stat[i].contam = 'clean'
    endif else begin
	stat[i].n_dark = 0
	stat[i].dark_mass = 0
        stat[i].contam = 'clean'
        stat[i].false = 'false'
    endelse	
  ENDFOR
ENDIF ELSE BEGIN
  dmass = double(d.mass[0])
  FOR i=0L, ngrps-1 do begin
    particles_considered = index[startind(amigaind[i]):endind(amigaind[i])]
    stat[i].gas_mass = 0
    stat[i].n_gas = 0
    stat[i].star_mass = 0
    stat[i].n_star = 0
    dark_ind = where(particles_considered LT h.ndark)
    stat[i].n_dark = n_elements(dark_ind) 
    stat[i].dark_mass = total(dmass[particles_considered(dark_ind)]) 
    test2 = where(dmass[particles_considered(dark_ind)] GT mindm)
    if test2[0] NE -1 then stat[i].contam = 'contam' else stat[i].contam = 'clean'
  ENDFOR
ENDELSE
"""


"""
; This next is designed to detect if the main halo has been broken down into smaller 
; false halos at the center, a common problem.
pot_false = where(stat.r_vir lt 5.0 and stat.n_star ne 0 and stat.n_dark lt 64, nfalse)
;pot_false = where(stat.r_vir lt 5.0 and stat.n_gas ne 0 and stat.n_star ne 0, nfalse)
;pot_false = where(stat.r_vir lt 1.0 and stat.n_dark lt 100, nfalse)
if nfalse ne 0 then begin 
    stat[pot_false].false = 'false' 
    stat[0].false = 'false?'
endif

;Find satellites
test_xc = stat.x_c-stat[0].x_c[0]
test_yc = stat.y_c-stat[0].y_c[0]
test_zc = stat.z_c-stat[0].z_c[0]
test_xc = test_xc[1:n_elements(test_xc)-1]
test_yc = test_yc[1:n_elements(test_xc)-1]
test_zc = test_zc[1:n_elements(test_xc)-1]
distance = ((test_xc^2.+test_yc^2.+test_zc^2.)^0.5)*1000. ;change from Mpc to kpc
sat = where(distance le stat[0].r_vir[0],comp=non)
stat[sat+1].sat = 'yes'
stat[non+1].sat = 'no'
stat[0].sat = 'no'
"""


"""
; Write the .stat file
get_lun, lun
filename = strcompress(tipsyfile[j]+'.amiga.stat')
openw, lun, filename
 printf, lun, format='(A5,2x,A9,2x,A8,2x,A8,2x,A8,2x,A16,2x,A10,2x,A16,2x,A16,2x,A16,2x,A10,2x,A9,2x,A10,2x,A9,2x,A9,2x,A9,2x,A10,2x,A10,2x,A10,2x,A7,2x,A10,2x,A7)', 'Grp', 'N_tot', 'N_gas', 'N_star', 'N_dark', 'Mvir(M_sol)', 'Rvir(kpc)', 'GasMass(M_sol)', 'StarMass(M_sol)', 'DarkMass(M_sol)', 'V_max', 'R@V_max', 'VelDisp', 'Xc', 'Yc', 'Zc', 'VXc', 'VYc', 'VZc', 'Contam', 'Satellite?', 'False?' ;'ID_A'

for i=0L,ngrps-1 do printf, lun, format='(I5,2x,I9,2x,I8,2x,I8,2x,I8,2x,A16,2x,D10.4,2x,A16,2x,A16,2x,A16,2x,D10.5,2x,D9.5,2x,D10.5,2x,D9.5,2x,D9.5,2x,D9.5,2x,D10.5,2x,D10.5,2x,D10.5,2x,A7,2x,A7,2x,A7)', stat[i].grp, stat[i].n_part, stat[i].n_gas, stat[i].n_star, stat[i].n_dark, stat[i].m_vir, stat[i].r_vir , stat[i].gas_mass*munit, stat[i].star_mass*munit, stat[i].dark_mass*munit, stat[i].v_max, stat[i].r_max, stat[i].sig_v, stat[i].x_c/1000., stat[i].y_c/1000., stat[i].z_c/1000., stat[i].vx_c, stat[i].vy_c, stat[i].vz_c, stat[i].contam, stat[i].sat, stat[i].false ;stat[i].amiga_orig_id

close, lun
free_lun, lun
"""


"""
; **************************************************************************
; Now, create the .gtp file
; **************************************************************************

statfile = strcompress(tipsyfile[j]+'.amiga.stat')
readcol, statfile, grp, npart, ngas, nstar, ndark, mvir, rvir, gmass, smass, dmass, vmax, rmax, sigv, xc, yc, zc, vxc, vyc, vzc, format='l,l,l,l,l,d,d,d,d,d,d,d,d,d,d,d,d,d,d'
  
star = replicate({mass:0.0, x:0.0, y:0.0, z:0.0, vx:0.0, vy:0.0, vz:0.0, metals:0.0, tform:0.0, eps:0.0, phi:0.0},n_elements(grp))
  for i=0L,n_elements(grp)-1 do begin
	star[i].mass = (gmass[i]+smass[i]+dmass[i])/munit
	star[i].x = (xc[i]*1000./(lunit))-0.5
	star[i].y = (yc[i]*1000./(lunit))-0.5
	star[i].z = (zc[i]*1000./(lunit))-0.5
	star[i].vx = vxc[i]/vunit
	star[i].vy = vyc[i]/vunit
	star[i].vz = vzc[i]/vunit
	star[i].eps = rvir[i]/lunit
  endfor
  star[*].metals = 0.0  
  star[*].phi = 0.0  
  star[*].tform = h.time
  header = h
  header.n = n_elements(grp)
  header.nstar = n_elements(grp)
  header.ngas = 0
  header.ndark = 0
  outfile = strcompress(tipsyfile[j]+'.amiga.gtp')
  wtipsy, outfile, header, gas, dark, star ;,/standard

ENDFOR

end
"""



