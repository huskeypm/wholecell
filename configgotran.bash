export MYPATH=/home/AD/pmke226/sources/mypython/
export LOC=/home/AD/pmke226/sources/
export PYTHONPATH=$PYTHONPATH:$MYPATH/lib/python2.7/site-packages/
export PATH=$PATH:$LOC/gotran/scripts/
python -c "import instant"
python -c "import modelparameters"
python -c "import gotran"

