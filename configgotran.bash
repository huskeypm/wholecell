
#export LOC=$HOME/sources/            
#export NAME=bdst227
export NAME=pmke226 
export LOC=/home/AD/$NAME/sources/            
export MYPATH=$LOC/mypython/
export PYTHONPATH=$PYTHONPATH:$MYPATH/lib/python2.7/site-packages/
export PATH=$PATH:$LOC/gotran/scripts/
python -c "import instant"
python -c "import modelparameters"
python -c "import gotran"

