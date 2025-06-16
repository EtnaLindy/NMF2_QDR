%% Nonnegative rank-2 factorization
% L*R' ~ N, computed with ANLS and 

% input 
% N:        a nonnegative matrix
% maxiter:  maximum number of iterations for ANLS (default 100)
% tol:      convergence criterion for ANLS (default 1e-3)

% output
% L:        mx2 nonnegative matrix
% R:        nx2 nonnegative matrix
% it:       ANLS iterations used
% conv:     convergence measure for consecutive ANLS iterates at the end


function [L,R,it,conv] = nmf2qdr(N,maxiter,tol)

if nargin < 2, maxiter = 100; end
if nargin < 3, tol = 1e-3; end

[L0,R0] = qdrinit(N);

[L,R,it,conv] = ANLS(N,L0,R0,maxiter,tol);


end
