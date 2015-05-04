#
#mpirun -np  8 python validationRoutines.py 1
#mpirun -np  8 python validationRoutines.py 2
#mpirun -np  8 python validationRoutines.py 3
#mpirun -np  8 python validationRoutines.py 4
#mpirun -np  8 python validationRoutines.py 4
#

import validationRoutines
validationRoutines.validation(1)
validationRoutines.validation(2)
validationRoutines.validation(3)
validationRoutines.validation(4)
validationRoutines.validation(5)
validationRoutines.validation(6)
validationRoutines.validation(7)             
validationRoutines.validation(13)
validationRoutines.validation(14)


print "All ValidationRoutines passed"
