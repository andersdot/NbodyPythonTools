import rockstarscript
rockstardir = '/home1/02575/lmanders/code/Rockstar-Galaxies/'
nnodes = 1
ncorespernode = 32
ServerInterface = 'ib0' #'ipogif0'
queue = 'largemem'

def make():
    rockstarscript.snaps()
    rockstarscript.cfg(ncorespernode=ncorespernode, nnodes=nnodes, ServerInterface=ServerInterface)
    rockstarscript.mainsubmissionscript(nnodes=nnodes, rockstardir=rockstardir, queue=queue)
    rockstarscript.postsubmissionscript(rockstardir=rockstardir, queue=queue)
    
#then submit <prefix to tipsyfile>.sbatch to the queue, and <prefixtotipsyfile>.post.sbatch to the queue depending on the prior finishing ok
#sbatch rockstar.sbatch
#sbatch --dependency=afterok:<jobid> rockstar.post.sbatch

if __name__=='__main__': make()
