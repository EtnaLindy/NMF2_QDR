
function [U,S,V] = fastsvd2(M)
    [m,n] = size(M);
    ratio = 10;

    if n >= ratio*m % thick/short
        MMt = M*M';
        [U,D] = eigs(MMt,2);
        S = sqrt(D);
        V = (M'*U)/S;
    elseif m >= ratio*n % tall/thin
        MtM = M'*M;
        [V,D] = eigs(MtM,2);
        S = sqrt(D);
        U = (M*V)/S;
    else
        [U,S,V] = svds(M,2);
    end

end