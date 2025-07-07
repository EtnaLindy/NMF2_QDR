%% SVD-based initialization for NMF2
% output L+ >= 0(mx2) and R+ >= 0(nx2) such that 
% L+*R+' = L*R', where L, R minimize
% ||Uhat - L||^2_2 + ||Vhat - R||_2^2
% where Uhat*Vhat' is the best rank-2 approximation for input N

% Assumes that N is nonnegative to begin with.

function [L,R] = qdrinit(N)

[U,S,V] = fastsvd2(N);

if S(2,2) == 0, L = U(:,1); R = V(:,1); return; end

U = U*sqrt(S); V = V*sqrt(S);

if any(U(:,1)<0), U=-U;V=-V; end % ensure that Perron vectors are nonnegative

% angular coordinates for U and V
psi=atan2(U(:,2),U(:,1)); du = sqrt(sum(U.^2,2));
phi=atan2(V(:,2),V(:,1)); dv = sqrt(sum(V.^2,2));

% Solve the problem by finding the zeros of the derivative piecewise
[alpha1,~] = mintheta(psi,phi,du,dv);
[alpha2,~] = mintheta(phi,psi,dv,du);

% find nonnegative L and R
s=sqrt(cos(alpha2-alpha1)); 
psi_new = min(max(psi,alpha1-pi/2),alpha2);
phi_new = min(max(phi,alpha2-pi/2),alpha1);
du = du.*cos(psi-psi_new);
dv = dv.*cos(phi-phi_new);
L=du.*[cos(alpha1-psi_new),sin(alpha2-psi_new)]/s;
R=dv.*[cos(alpha2-phi_new),sin(alpha1-phi_new)]/s;
D=diag(diag(L'*L)./diag(R'*R))^(1/4); %balance
L=L/D;R=R*D;


end

