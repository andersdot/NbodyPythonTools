import nbdpt.rockstar.rockstarscript as rs

rockstardir   = '/home1/02575/lmanders/code/Rockstar-Galaxies/'
nnodes        = 1
ncorespernode = 32 
queue         = 'largemem'                 #largemem, normal on stampede
email         = 'l.sonofanders@gmail.com'  #please change so I don't get all your emails :P
machine       = 'stampede'                 #stampede or pleiades
walltime      = '24:00:00'                 #need to use this notation 
massdef       = '200c'                     #mass options, 'vir', '###b', '###c'
massdef2      = None
ServerInterface = 'ib0'                    #'ipogif0' on bluewaters

def make():
    rs.snaps()
    rs.cfg(ncorespernode=ncorespernode, nnodes=nnodes, 
           ServerInterface=ServerInterface, massdef=massdef, 
           massdef2=massdef2) 
    rs.mainsubmissionscript(nnodes=nnodes, ncorespernode=ncorespernode, 
                            machine=machine, email=email, 
                            rockstardir=rockstardir, queue=queue)
    rs.postsubmissionscript(nnodes=nnodes, ncorespernode=ncorespernode,
                            machine=machine, email=email, 
                            rockstardir=rockstardir, queue=queue)
    
#then submit rockstar.sbatch to the queue, and rockstar.post.sbatch to the queue 
#    depending on the prior finishing ok

#sbatch rockstar.sbatch
#sbatch --dependency=afterok:<jobid> rockstar.post.sbatch

if __name__=='__main__': make()
