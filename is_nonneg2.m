function [is_nonneg] = is_nonneg2(N)
    
    [U,S,V] = svds(N,2);
    if any(U(:,1) < 0), U(:,1) = -U(:,1); V(:,1) = -V(:,1); end
    
     % we hope that Uhat(:,1) does not have zero elements ehehe
     s = sqrt(S(2,2)/S(1,1));
     Urat = s*U(:,2)./U(:,1);
     Vrat = s*V(:,2)./V(:,1);
     is_nonneg =  min(Urat)*max(Vrat) >= -1 && max(Urat)*min(Vrat) >= -1;

end