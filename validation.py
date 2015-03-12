#
#mpirun -np  8 python validations.py 1
#mpirun -np  8 python validations.py 2
#mpirun -np  8 python validations.py 3
#mpirun -np  8 python validations.py 4
#mpirun -np  8 python validations.py 4
#

import validations
validations.validation(1)
validations.validation(2)
validations.validation(3)
validations.validation(4)
validations.validation(5)


print "All Validations passed"
