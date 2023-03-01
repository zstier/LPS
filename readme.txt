Python3 implementations of LPS Ramanujan graphs and associated walk operators
Zachary Stier
updated February 2023

These Python files modularize the construction of the simple random walk (SRW) and non-backtracking random walk (NBRW) operators for LPS Ramanujan graphs, as well as the SRW operator for the LPS quotient graphs Y. Most files work by pickle-ing a specific computation result, to be called by a subsequent file. As such there is a specific workflow. Running "python3 setupLPS.py p q" does this workflow all at once. 

It is important that the directories jacobi4, srwadjmat, srwadjlist, and nbrwadjmat are all present. 

"python3 findJacobi4.py n" finds all quadruples (a,b,c,d) whose sum of squares is n, with a positive and odd. pickle-d output is "jacboi4/n.p".

"python3 makeLPSadjlist.py p q" creates an adjacency list (dictionary) representing the LPS graph X^p,q. Requires "jacobi4/p.p". pickle-d output is "srwadjlist/p,q.p". 

"python3 makeLPSadjmat.py p q" creates an adjacency matrix (numpy array) representing the LPS graph X^p,q. Requires "jacobi4/p.p". pickle-d output is "srwadjmat/p,q.p". 

"python3 makeNBRWfromLPS.py p q" creates a csr_matrix representing the nonnormalized NBRW on the LPS graph X^p,q (so that each row sum is exactly p). Requires "srwadjlist/p,q.p". pickle-d output is "nbrwadjmat/p,q.p". 

"python3 makeYadjmat.py p q" creates an adjacency matrix (array) representing the LPS graph Y^p,q. Requires "jacobi4/p.p". pickle-d output is "Y/p,q.p". 

"python3 setup_small_sparse_Ys.py" creates the adjacency matrix (array) for each of 97 small LPS Y graphs. 