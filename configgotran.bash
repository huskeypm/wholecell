export MYPATH=/home/huskeypm/sources/mypython/
export LOC=/home/huskeypm/sources/
export PYTHONPATH=$PYTHONPATH:$MYPATH/lib/python2.7/site-packages/

python -c "import instant"
python -c "import modelparameters"
python -c "import gotran"
export PYTHONPATH=$PYTHONPATH:/home/huskeypm/sources/mypython//lib/python2.7/site-packages/
export PATH=$PATH:$LOC/gotran/scripts/

