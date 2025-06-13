function [L,R,it,conv] = nmf2qdr(N,maxiter,tol)

if nargin < 2, maxiter = 100; end
if nargin < 3, tol = 1e-3; end

[L0,R0] = qdrinit(N);

[L,R,it,conv] = ANLS(N,L0,R0,maxiter,tol);


end
