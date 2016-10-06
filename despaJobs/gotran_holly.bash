#!/bin/bash
#PBS -l nodes=1:ppn=64 
#PBS -q long    
#PBS -m abe

## For running dolfin jobs on holly 
### NOTE: Run 'module load FEniCS' interactively before
### using this with qsub.
export BASEDIR=$PBS_O_WORKDIR
export PROCS=64
export OMPI_MCA_orte_default_hostfile=$PBS_NODEFILE
export OMPI_MCA_orte_leave_session_attached=1
. /etc/profile.d/modules.sh

module load FEniCS.15
export LOC=$HOME/sources
export MYPATH=$LOC/mypython
export PYTHONPATH=$PYTHONPATH:$MYPATH/lib/python2.7/site-packages/


cd /home/pmke226/sources/wholecell/



python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 1.000000 -name ./run_G_CaBk1.00_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 1.250000 -name ./run_G_CaBk1.00_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 1.500000 -name ./run_G_CaBk1.00_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 1.750000 -name ./run_G_CaBk1.00_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 2.000000 -name ./run_G_CaBk1.00_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 2.250000 -name ./run_G_CaBk1.00_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 2.500000 -name ./run_G_CaBk1.00_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.000000 -var G_NaBk 2.750000 -name ./run_G_CaBk1.00_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 1.000000 -name ./run_G_CaBk1.25_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 1.250000 -name ./run_G_CaBk1.25_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 1.500000 -name ./run_G_CaBk1.25_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 1.750000 -name ./run_G_CaBk1.25_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 2.000000 -name ./run_G_CaBk1.25_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 2.250000 -name ./run_G_CaBk1.25_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 2.500000 -name ./run_G_CaBk1.25_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.250000 -var G_NaBk 2.750000 -name ./run_G_CaBk1.25_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 1.000000 -name ./run_G_CaBk1.50_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 1.250000 -name ./run_G_CaBk1.50_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 1.500000 -name ./run_G_CaBk1.50_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 1.750000 -name ./run_G_CaBk1.50_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 2.000000 -name ./run_G_CaBk1.50_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 2.250000 -name ./run_G_CaBk1.50_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 2.500000 -name ./run_G_CaBk1.50_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.500000 -var G_NaBk 2.750000 -name ./run_G_CaBk1.50_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 1.000000 -name ./run_G_CaBk1.75_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 1.250000 -name ./run_G_CaBk1.75_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 1.500000 -name ./run_G_CaBk1.75_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 1.750000 -name ./run_G_CaBk1.75_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 2.000000 -name ./run_G_CaBk1.75_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 2.250000 -name ./run_G_CaBk1.75_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 2.500000 -name ./run_G_CaBk1.75_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 1.750000 -var G_NaBk 2.750000 -name ./run_G_CaBk1.75_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 1.000000 -name ./run_G_CaBk2.00_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 1.250000 -name ./run_G_CaBk2.00_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 1.500000 -name ./run_G_CaBk2.00_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 1.750000 -name ./run_G_CaBk2.00_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 2.000000 -name ./run_G_CaBk2.00_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 2.250000 -name ./run_G_CaBk2.00_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 2.500000 -name ./run_G_CaBk2.00_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.000000 -var G_NaBk 2.750000 -name ./run_G_CaBk2.00_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 1.000000 -name ./run_G_CaBk2.25_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 1.250000 -name ./run_G_CaBk2.25_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 1.500000 -name ./run_G_CaBk2.25_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 1.750000 -name ./run_G_CaBk2.25_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 2.000000 -name ./run_G_CaBk2.25_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 2.250000 -name ./run_G_CaBk2.25_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 2.500000 -name ./run_G_CaBk2.25_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.250000 -var G_NaBk 2.750000 -name ./run_G_CaBk2.25_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 1.000000 -name ./run_G_CaBk2.50_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 1.250000 -name ./run_G_CaBk2.50_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 1.500000 -name ./run_G_CaBk2.50_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 1.750000 -name ./run_G_CaBk2.50_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 2.000000 -name ./run_G_CaBk2.50_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 2.250000 -name ./run_G_CaBk2.50_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 2.500000 -name ./run_G_CaBk2.50_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.500000 -var G_NaBk 2.750000 -name ./run_G_CaBk2.50_G_NaBk2.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 1.000000 -name ./run_G_CaBk2.75_G_NaBk1.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 1.250000 -name ./run_G_CaBk2.75_G_NaBk1.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 1.500000 -name ./run_G_CaBk2.75_G_NaBk1.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 1.750000 -name ./run_G_CaBk2.75_G_NaBk1.75_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 2.000000 -name ./run_G_CaBk2.75_G_NaBk2.00_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 2.250000 -name ./run_G_CaBk2.75_G_NaBk2.25_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 2.500000 -name ./run_G_CaBk2.75_G_NaBk2.50_stim1000 &
 python runShannonTest.py -stim 1000 -T 300000 -var G_CaBk 2.750000 -var G_NaBk 2.750000 -name ./run_G_CaBk2.75_G_NaBk2.75_stim1000 &
sleep 150000 
print "LEAVING!" 

