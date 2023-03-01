import os
from os.path import exists
import sys

p, q = 5, 29

if len(sys.argv) > 2:
	p = int(sys.argv[1])
	q = int(sys.argv[2])

if not exists("jacobi4/" + str(p) + ".p"):
	os.system("python3 findJacobi4.py " + str(p))
	print("found Jacobi quadruples")
if not exists("srwadjlist/" + str(p) + "," + str(q) + ".p"):
	os.system("python3 makeLPSadjlist.py " + str(p) + " " + str(q))
	print("created SRW adjacency list")
if not exists("srwadjmat/" + str(p) + "," + str(q) + ".p"):
	os.system("python3 makeLPSadjmat.py " + str(p) + " " + str(q))
	print("created SRW adjacency matrix")
if (not exists("nbrwadjmat/" + str(p) + "," + str(q) + ".p")) or (not exists("nbrwadjlist/" + str(p) + "," + str(q) + ".p")):
	os.system("python3 makeNBRWfromLPS.py " + str(p) + " " + str(q))
	print("created NBRW adjacency matrix/list")
