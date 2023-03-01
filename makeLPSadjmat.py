import pickle
import numpy as np
import sys

def legendre(a, p): # Legendre symbol
	return pow(a, (p - 1) // 2, p)

def tonelli(n, p): # Tonelli-Shanks algorithm, from https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python
	assert legendre(n, p) == 1, "not a square mod p"
	q = p - 1
	s = 0
	while q % 2 == 0:
		q //= 2
		s += 1
	if s == 1:
		return pow(n, (p + 1) // 4, p)
	for z in range(2, p):
		if p - 1 == legendre(z, p):
			break
	c = pow(z, q, p)
	r = pow(n, (q + 1) // 2, p)
	t = pow(n, q, p)
	m = s
	t2 = 0
	while (t - 1) % p != 0:
		t2 = (t * t) % p
		for i in range(1, m):
			if (t2 - 1) % p == 0:
				break
			t2 = (t2 * t2) % p
		b = pow(c, 1 << (m - i - 1), p)
		r = (r * b) % p
		c = (b * b) % p
		t = (t * c) % p
		m = i
	return r

p, q = 5, 29
iota = 12 # sqrt(-1) mod q
sqrtpinv = 8 # sqrt(p) mod q

if len(sys.argv) > 2:
	p = int(sys.argv[1])
	q = int(sys.argv[2])
	iota = tonelli(-1, q)
	sqrtpinv = tonelli(pow(p,-1,q), q)

pquadruples = pickle.load(open("jacobi4/" + str(p)+".p", "rb")) # all quadruples (a,b,c,d) with a positive and odd and sum of squares equalling p

A = [((sqrtpinv * (a+b*iota)) % q, (sqrtpinv * (c+d*iota)) % q, (sqrtpinv * (-c+d*iota)) % q, (sqrtpinv * (a-b*iota)) % q) for (a,b,c,d) in pquadruples] # identify with elements of PSL2(Fq)

V = [] # vertex set of X^{p,q}
### "standard form" for PSL2(Fq): flatten ((a,b),(c,d)) where a < q//2 and if a = 0 then b < q//2
for b in range(1,1+q//2): # a = 0
	c = pow(-b,-1,q)
	for d in range(q):
		V.append((0,b,c,d))
for a in range(1,1+q//2): # a > 0
	for b in range(q):
		for c in range(q):
			d = (1+b*c)*pow(a,-1,q) % q
			V.append((a,b,c,d))

# print(len(V))

def qmul(X,Y): # multiply tuples (flattened 2x2 mats) mod q and use the standard form
	(a,b,c,d),(e,f,g,h) = X,Y
	(i,j,k,l) = ((a*e+b*g) % q, (a*f+b*h) % q, (c*e+d*g) % q, (c*f+d*h) % q)
	if i == 0:
		if j > q//2:
			return (0, -j % q, -k % q, -l % q)
		return (0, j, k, l)
	if i > q//2:
		return (-i % q, -j % q, -k % q, -l % q)
	return (i, j, k, l)

def n(X): # "integer form" of a tuple
### this method may be a problem but for now it seems to work fine
	(a,b,c,d) = X
	return a*q**3 + b*q**2 + c*q + d

N = len(V)
hit = [False]*(q**4)
for v in V:
	hit[n(v)] = True
ind = [0]*(q**4)
tot = 0
for i in range(q**4):
	ind[i] = tot
	if hit[i]:
		tot += 1

adj = np.zeros((N,N)) # adjacency matrix; work with n(v) of matrices
for v in V:
	nv = n(v)
	for a in A:
		va = qmul(v,a)
		nva = n(va)
		adj[ind[nv]][ind[nva]] = 1.0

pickle.dump(adj, open("srwadjmat/" + str(p) + "," + str(q) + ".p", "wb"))

