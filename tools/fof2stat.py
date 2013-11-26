import nptipsyreader
import glob
import pdb
import numpy as np

def fofgrp2stat(grpfilename, gtpfilename, writestat=True):
    grpf = open(grpfilename)
    statfilename = ('.').join(grpfilename.split('.')[:-1]) + '.stat'
    statdtype = np.dtype({'names': ('grp','mvir[Msol]','npart','ngas','nstar','ndark',
                                    'gas[Msol]','star[Msol]','dark[Msol]', 'rmax[kpc]',
                                    'x[Mpc]', 'y[Mpc]', 'z[Mpc]', 'vx[km/s]', 'vy[km/s]', 'vz[km/s]'),
                          'formats': ('int32'  ,'float64','int32'  ,'int32'  ,'int32'  ,'int32'  ,
                                      'float64','float64','float64','float64','float64','float64',
                                      'float64','float64','float64','float64','float64','float64',
                                      'float64','int32'  ,'bool'   ,'bool'   ,'bool')})
    """
    statdtype = np.dtype({'names': ('grp','mvir[Msol]','npart','ngas','nstar','ndark',
                                    'gasmass[Msol]','starmass[Msol]','darkmass[Msol]',
                                    'rvir[kpc]', 'vmax[km/s]', 'rmax[kpc]', 'sigv[km/s]',
                                    'x[Mpc]', 'y[Mpc]', 'z[Mpc]', 'vx[km/s]', 'vy[km/s]', 'vz[km/s]', 
                                    'amorigid','sat', 'contam', 'false'),
                          'formats': ('int32'  ,'float64','int32'  ,'int32'  ,'int32'  ,'int32'  ,
                                      'float64','float64','float64','float64','float64','float64',
                                      'float64','float64','float64','float64','float64','float64',
                                      'float64','int32'  ,'bool'   ,'bool'   ,'bool')})

    """
    grp = [line.split('\n')[0] for line in grpf]
    grp = np.array(grp, dtype='int32')
    grp = grp[1:]
    grpf.close()

    gtp = nptipsyreader.Tipsy(gtpfilename)
    gtp._read()

    tipsy = nptipsyreader.Tipsy(('.').join(grpfilename.split('.')[:-4]))
    tipsy._read_param()
    tipsy._read()

    massunit = np.float(tipsy.paramfile['dMsolUnit'])
    kpcunit  = np.float(tipsy.paramfile['dKpcUnit'])
    velunit  = tipsy.velunit

    halos = np.unique(grp)[1:]
    stat = np.empty(len(halos), dtype=statdtype)

    if writestat: 
        statf = open(statfilename, 'w')
        statf.write('# ' + "{0[0]:>10}  {0[1]:>10}  {0[2]:>10}  {0[3]:>10}  {0[4]:>10}  {0[5]:>10}  {0[6]:>10}  {0[7]:>10}  {0[8]:>10}  {0[9]:>10}  {0[10]:>10}  {0[11]:>10}  {0[12]:>10}  {0[13]:>10}  {0[14]:>10}  {0[15]:>10}".format(stat.dtype.names) + '\n')
                                                                            
    for i in range(len(halos)):
        thishalo = (grp == halos[i])
        thishalogas  = thishalo[0                      : tipsy.ngas-1                        ]
        thishalodark = thishalo[tipsy.ngas             : tipsy.ngas+tipsy.ndark-1            ]
        thishalostar = thishalo[tipsy.ngas+tipsy.ndark : tipsy.ngas+tipsy.ndark+tipsy.nstar-1]

        stat[i]['grp']    = halos[i]
        #stat[i]['sat']    = False
        #stat[i]['contam'] = False
        #stat[i]['false']  = False
        stat[i]['npart']  = np.sum(thishalo).astype('int32')

        stat[i]['ngas']  = np.sum(thishalogas ).astype('int32')
        stat[i]['ndark'] = np.sum(thishalodark).astype('int32')
        stat[i]['nstar'] = np.sum(thishalostar).astype('int32')

        stat[i]['mvir[Msol]']  = (np.sum(tipsy.gas ['mass'][thishalogas] ).astype('float64') + 
                            np.sum(tipsy.dark['mass'][thishalodark]).astype('float64') + 
                            np.sum(tipsy.star['mass'][thishalostar]).astype('float64'))*massunit
        
        #stat[i]['rvir[kpc]'] = False
        
        stat[i]['gas[Msol]']  = np.sum(tipsy.gas['mass'] [thishalogas] ).astype('float64')*massunit
        stat[i]['dark[Msol]'] = np.sum(tipsy.dark['mass'][thishalodark]).astype('float64')*massunit
        stat[i]['star[Msol]'] = np.sum(tipsy.star['mass'][thishalostar]).astype('float64')*massunit

        #stat[i]['vmax[km/s]'] = False
        #stat[i]['sigv[km/s]'] = False
        stat[i]['x[Mpc]'] = gtp.star[i]['x']*np.float(tipsy.paramfile['dKpcUnit'])/1000.
        stat[i]['y[Mpc]'] = gtp.star[i]['y']*np.float(tipsy.paramfile['dKpcUnit'])/1000.
        stat[i]['z[Mpc]'] = gtp.star[i]['z']*np.float(tipsy.paramfile['dKpcUnit'])/1000.
        stat[i]['vx[km/s]'] = gtp.star[i]['vx']*tipsy.velunit
        stat[i]['vy[km/s]'] = gtp.star[i]['vy']*tipsy.velunit
        stat[i]['vz[km/s]'] = gtp.star[i]['vz']*tipsy.velunit
        stat[i]['rmax[kpc]'] = gtp.star[i]['eps']*np.float(tipsy.paramfile['dKpcUnit'])
        
        #stat[i]['amorigid'] = False
        #import IPython; IPython.embed()

    stat = stat[np.argsort(stat['mvir[Msol]'])[::-1]]
    for i in range(len(stat)):

        if writestat: statf.write('  ' + "{stat[0]:>10n}  {stat[1]:>10.4e}  {stat[2]:>10n}  {stat[3]:>10n}  {stat[4]:>10n}  {stat[5]:>10n}  {stat[6]:>10.4e}  {stat[7]:>10.4e}  {stat[8]:>10.4e}  {stat[9]:>10.2f}  {stat[10]:>10.2f}  {stat[11]:>10.2f}  {stat[12]:>10.2f}  {stat[13]:>10.2f}  {stat[14]:>10.2f}  {stat[15]:>10.2f}".format(stat=stat[i]) + '\n')


    statf.close()
    return stat

def allfof2stat():
    grpfilenames = glob.glob('*.grp') #'*.fof*.grp')
    gtpfilenames = glob.glob('*.gtp')#('*.fof*.gtp')[0]
    grpfilenames.sort()
    gtpfilenames.sort()

    alreadymade_stat = glob.glob('*.fof.*.stat')
    alreadymade_stat.sort()

    found       = [('.').join(filename.split('.')[:-4]) for filename in grpfilenames]
    alreadymade = [('.').join(filename.split('.')[:-4]) for filename in alreadymade_stat]

    wanted = list(set(alreadymade) ^ set(found))
    wanted.sort()

    for i in range(len(wanted)):
        grpfilename = glob.glob(wanted[i] + '.fof.*.grp')[0]
        gtpfilename = glob.glob(wanted[i] + '.fof.*.gtp')[0]
        print grpfilename, gtpfilename
        fofgrp2stat(grpfilename, gtpfilename, writestat=True)

