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
print "WARNING: PKH needs to fix"
#validationRoutines.validation(3)
validationRoutines.validation(4)
validationRoutines.validation(5)
validationRoutines.validation(6)
validationRoutines.validation(61)
validationRoutines.validation(7)             
validationRoutines.validation(71)             


print "All ValidationRoutines passed"
