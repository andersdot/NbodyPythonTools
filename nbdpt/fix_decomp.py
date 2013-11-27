
def fix_decomp():

    """
    fix_decomp, 'cosmo6.25PLK.192g1bwK1C52.003697.amiga.grp', 'cosmo6.25PLK.192g1bwK1C52.003697.amiga.grp', 'cosmo6.25PLK.192g1bwK1C52.003697.h1.decomp', 'cosmo6.25PLK.192g1bwK1C52.003697.h1.decomp.fixed'
    
    grp= read_lon_array(grpfile)
    rtipsy, tipsyfile, h,g,d,s
    totdg = h.ngas+h.ndark
    ind=lindgen(h.n)
    w=where(grp eq 190 and ind ge totdg)
    ;read in the new stupid decomp file (boo on the brain for this!)
    rdfloat, decompfile, mark
    ngrp=n_elements(grp)
    new=lonarr(ngrp)
    stop
    new[w]=mark
    openw, lun, outfile, /get_lun
    printf, lun, ngrp
    for j=0L, ngrp-1 do printf, lun, new[j]
    close, lun
    free_lun, lun
    stop
    end
    """
