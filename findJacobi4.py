import pickle
import numpy as np
from numpy import sqrt
import sys

p = 5

if len(sys.argv) > 1:
	p = int(sys.argv[1])

A = [] # list of quadruples

for a in range(1, int(sqrt(p))+1, 2): # look for odd a
	rb = int(sqrt(p-a**2))
	fb = rb % 2
	for b in range(-rb+fb, rb+1, 2): # look for even b which are not too large
		rc = int(sqrt(p-a**2-b**2))
		fc = rc % 2
		for c in range(-rc+fc, rc+1, 2): # look for even c which are not too large
			d = int(sqrt(p-a**2-b**2-c**2))
			if a**2 + b**2 + c**2 + d**2 == p: # see if (a,b,c,d) is viable
				A.append((a,b,c,d))
				if d != 0:
					A.append((a,b,c,-d))

pickle.dump(A, open("jacobi4/" + str(p) + ".p", "wb"))
