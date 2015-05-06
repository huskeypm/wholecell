#!/bin/bash
#PBS -l nodes=1:ppn=10 
#PBS -q long    
#PBS -m abe

## For running dolfin jobs on holly 
### NOTE: Run 'module load FEniCS' interactively before
### using this with qsub.
export BASEDIR=$PBS_O_WORKDIR
export PROCS=32
export OMPI_MCA_orte_default_hostfile=$PBS_NODEFILE
export OMPI_MCA_orte_leave_session_attached=1
. /etc/profile.d/modules.sh

module load FEniCS.15
export LOC=$HOME/sources
export MYPATH=$LOC/mypython
export PYTHONPATH=$PYTHONPATH:$MYPATH/lib/python2.7/site-packages/


cd /home/pmke226/sources/wholecell/

python runShannonTest.py -stim 1000 -T 300000 -name /home/pmke226/scratch/despa/despahealthy1000_300.pickle & 
python runShannonTest.py -stim 1000 -T 250000 -name /home/pmke226/scratch/despa/despahealthy1000_250.pickle & 
python runShannonTest.py -stim 1000 -T 200000 -name /home/pmke226/scratch/despa/despahealthy1000_200.pickle & 
python runShannonTest.py -stim 1000 -T 150000 -name /home/pmke226/scratch/despa/despahealthy1000_150.pickle & 
python runShannonTest.py -stim 1000 -T 100000 -name /home/pmke226/scratch/despa/despahealthy1000_100.pickle & 
python runShannonTest.py -stim 1000 -T 75000 -name /home/pmke226/scratch/despa/despahealthy1000_75.pickle & 
python runShannonTest.py -stim 1000 -T 50000 -name /home/pmke226/scratch/despa/despahealthy1000_50.pickle & 
python runShannonTest.py -stim 1000 -T 25000 -name /home/pmke226/scratch/despa/despahealthy1000_25.pickle & 
python runShannonTest.py -stim 1000 -T 15000 -name /home/pmke226/scratch/despa/despahealthy1000_15.pickle & 
python runShannonTest.py -stim 1000 -T 10000 -name /home/pmke226/scratch/despa/despahealthy1000_10.pickle & 


sleep 150000 
print "LEAVING!" 

