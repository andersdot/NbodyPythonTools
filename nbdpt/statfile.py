import numpy as np
import pdb

#function read_stat_struc_AMIGA, file
#FMT = 'F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,X'
#readcol, file, F=FMT, group, member, ngas, nstar, n_dark, m_tot, R_vir, $
#  m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
#  xVcm, yVcm, zVcm, /silent
#halos = replicate({group:0, npart:0L, m_tot:0., m_gas:0., m_star:0., $
#                  vc:0., r_vc_max:0., ngas:0L, nstar:0L, rvir:0., $
#                  v_disp:0., xc:0., yc:0., zc:0., xVcm:0., $
#                   yVcm:0., zVcm:0.}, n_elements(group))
#dtype = 'F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,X,X'


def read(file):
	dtype={'names': ('group', 'ntot', 'ngas', 'nstar', 'ndark', 'mtot', 
			 'Rvir', 'mgas', 'mstar', 'mdark', 'vcmax', 'rvcmax', 
			 'vdisp',  'xc', 'yc', 'zc', 'xVcm', 'yVcm', 'zVcm', 'Contam', 'satellite', 'false'),
	       'formats': ('i4','i4','i4','i4','i4','f4',
			   'f4','f4','f4','f4','f4','f4',
			   'f4','f4','f4','f4','f4','f4','f4', 'S6', 'S6', 'S6')}
	data = np.loadtxt(file, dtype=dtype, skiprows=1)

	return data

def readfof(file):
	dtype = {'names': ('grp','mtot','npart','ngas','nstar','ndark',
			   'gas','star','dark','rmax','x','y','z','vx','vy','vz'),
		 'formats':('i4','f4','i4', 'i4', 'i4', 'i4',
			    'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
	

	data = np.loadtxt(file, dtype=dtype, skiprows=1)

	return data
#readcol, file, F=FMT, group, ntot, ngas, nstar, n_dark, m_tot, R_vir, $
#  m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
#  xVcm, yVcm, zVcm, /silent

#halos = replicate({group:0, npart:0, m_tot:0., m_gas:0., m_star:0., $
#                  vc:0., r_vc_max:0., ngas:0L, nstar:0L, rvir:0., $
#                  v_disp:0., xc:0., yc:0., zc:0., xVcm:0., $
#                   yVcm:0., zVcm:0.}, n_elements(group))
"""
halos.vc = REFORM(vc_max)
halos.group = REFORM(group)
halos.npart = REFORM(ntot)
halos.ngas = REFORM(ngas)
halos.nstar = REFORM(nstar)
halos.m_tot = REFORM(m_tot)
halos.m_gas = REFORM(m_gas)
halos.m_star = REFORM(m_star)
halos.rvir = REFORM(R_vir)
halos.r_vc_max = REFORM(r_vc_max)
halos.v_disp = REFORM(v_disp)
halos.xc = REFORM(xc)
halos.xVcm = REFORM(xVcm)
halos.yc = REFORM(yc)
halos.yVcm = REFORM(yVcm)
halos.zc = REFORM(zc)
halos.zVcm = REFORM(zVcm)

RETURN, halos
END
"""
