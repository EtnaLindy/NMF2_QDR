# NMF2_QDR

MATLAB algorithm for computing a nonnegative rank-2 approximation for the given nonnegative matrix.

This can be done using the command nmg2qdr(N) where N is a mxn nonnegative matrix. 
The output is L (mx2) and R (nx2) such that L*R' gives a nonnegative rank-2 approximation for N.

The main contribution is the initialization method qdrinit(N).


Includes benchmarking the initialization method against other methods for computing NMF2 approximations. 
The raw data is in RESULTS folder, and the files benchmarking_i.m for i in {a,b,c,d} can be used to produce 
figures for the various test cases.

To rerun the benchmarking computations, one needs the codes for SPA, NNSVDLRC, and NNDSVD. All of these can 
be found in https://sites.google.com/site/nicolasgillis/code.

