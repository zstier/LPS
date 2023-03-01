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
# goal: create the LPS graph Y^p,q

p, q = 17, 13
iota = 12 # sqrt(-1) mod q

if len(sys.argv) > 2:
	p = int(sys.argv[1])
	q = int(sys.argv[2])
	iota = tonelli(-1, q)

V = range(q+1)

A = [[0 for _ in range(q+1)] for _ in range(q+1)]
pquadruples = pickle.load(open("jacobi4/" + str(p)+".p", "rb")) # all quadruples (a,b,c,d) with a positive and odd and sum of squares equalling p

for v in V:
	for (a,b,c,d) in pquadruples:
		w = (a+b*iota) % q
		x = (c+d*iota) % q
		y = (-c+d*iota) % q
		z = (a-b*iota) % q
		if v == q:
			num = w
			den = y
		else:
			num = (w*v + x) % q
			den = (y*v + z) % q
		if den == 0:
			A[v][q] += 1
		else:
			t = (num*pow(den,-1,q)) % q
			A[v][t] += 1
			
pickle.dump(A, open("Y/" + str(p) + "," + str(q) + ".p", "wb"))
