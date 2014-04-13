import nbdpt.rockstar.rockstarscript as rs

rockstardir   = '/home1/02575/lmanders/code/Rockstar-Galaxies/'
nnodes        = 1
ncorespernode = 32 
queue         = 'largemem'                 #largemem, normal on stampede
email         = 'l.sonofanders@gmail.com'  #please change so I don't get all your emails :P
machine       = 'stampede'                 #stampede, pleiades, bluewaters or interactive
walltimemain  = '24:00:00'                 #need to use this notation 
walltimepost  = '3:00:00'
massdef       = '200c'                     #mass options, 'vir', '###b', '###c'
massdef2      = None
ServerInterface = 'ib0'                    #'ipogif0' on bluewaters
fileformat = 'TIPSY'

def make():
    rs.snaps()
    rs.cfg(ncorespernode=ncorespernode, nnodes=nnodes, 
           ServerInterface=ServerInterface, massdef=massdef, 
           massdef2=massdef2, fileformat=fileformat) 
    rs.mainsubmissionscript(nnodes=nnodes, ncorespernode=ncorespernode, 
                            machine=machine, email=email, 
                            rockstardir=rockstardir, queue=queue, walltime=walltimemain)
    rs.postsubmissionscript(nnodes=nnodes, ncorespernode=ncorespernode,
                            machine=machine, email=email, 
                            rockstardir=rockstardir, queue=queue, walltime=walltimepost, fileformat=fileformat)
    
#then submit rockstar.sbatch to the queue, and rockstar.post.sbatch to the queue 
#    depending on the prior finishing ok

#sbatch rockstar.sbatch
#sbatch --dependency=afterok:<jobid> rockstar.post.sbatch

if __name__=='__main__': make()
