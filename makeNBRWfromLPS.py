import pickle
import numpy as np
import sys
from scipy.sparse import csr_matrix

p, q = 5, 29

if len(sys.argv) > 2:
	p = int(sys.argv[1])
	q = int(sys.argv[2])

al = pickle.load(open("srwadjlist/" + str(p) + "," + str(q) + ".p", "rb")) # adjacency list

d = p+1 # degree
n = (q-1)*q*(q+1)//2 # number of vertices
e = n*d # number of arcs

verts = sorted(al.keys())
inverse = {} # if verts[i] = j then inverse[j] = i
for i in range(n):
	j = verts[i]
	inverse[j] = i

D = {} # initialize dict for adjacency list encoding
row_ind = [] # initialize list for sparse matrix encoding
col_ind = [] # initialize list for sparse matrix encoding
for v in range(n):
	for i in range(d):
		D[d*v+i] = []
		x = verts[v] # current vertex, identified as a base-q matrix (its inverse is i)
		y = al[x][i] # ith neighbor of x, identified as a base-q matrix
		if y not in al.keys(): # hopefully should never occur!
			print("uh oh!")
			print(al[x])
			print(x,y)
		w = inverse[y] # y's sequential position
		for j in range(d):
			z = al[y][j] # jth neighbor of y
			if z != x: # can't backtrack
				row_ind.append(d*v+i)
				col_ind.append(d*w+j)
				D[d*v+i].append(d*w+j)

data = [1.0]*(e*(d-1)) # just put 1's in each entry

B = csr_matrix((data, (row_ind, col_ind))) # create sparse NBRW matrix

pickle.dump(B, open("nbrwadjmat/" + str(p) + "," + str(q) + ".p", "wb"))
pickle.dump(D, open("nbrwadjlist/" + str(p) + "," + str(q) + ".p", "wb"))
pickle.dump(inverse, open("inverse/" + str(p) + "," + str(q) + ".p", "wb"))
