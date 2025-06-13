function [L,R] = qdrinit(N)
%QDRINIT Summary of this function goes hereDu

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

s=sqrt(cos(alpha2-alpha1)); 
L=du.*[cos(alpha1-psi),sin(alpha2-psi)]/s;
R=dv.*[cos(alpha2-phi),sin(alpha1-phi)]/s;
D=diag(diag(L'*L)./diag(R'*R))^(1/4);
L=L/D;R=R*D;

end

