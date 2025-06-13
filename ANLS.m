function [L,R,it,conv] = ANLS(M,L0,R0,maxiter,tol)

if nargin < 4, maxiter = 100; end
if nargin < 5, tol = 1e-3; end 

L = L0; R = R0;
Mt = M';

j = 0;
conv = 1;

while j < maxiter && conv > tol
    R_new = NNLS(M,L);
    conv = sum(sqrt(sum((R-R_new).^2,2)./sqrt(sum(R_new.^2,2))));
    R = R_new;

    L_new = NNLS(Mt,R);
    conv = conv + sum(sqrt(sum((L-L_new).^2,2)./sqrt(sum(L_new.^2,2)))); 
    L = L_new;
    j = j + 1;
end

if conv > tol, j = j + 1; end

it = j;


end

