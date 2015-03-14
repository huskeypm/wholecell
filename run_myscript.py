#!/usr/bin/env python

#import numpy as np
import os

#iters = 3
#phis=np.linspace(0.1,1.0,iters)
#Kds=np.linspace(-6,-4,iters)

phis=0.1,0.55,1.0
Kds=-6,-5,-4

for i,Kdi in enumerate(Kds):
    for j,phij in enumerate(phis):
	os.system("mpirun -np 30 python myscript.py -runMPI asdf -Kdi Kdi -phij phij")
