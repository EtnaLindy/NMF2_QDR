function M2 = eym(U,r)

[L,S,R] = svds(U,r);
M2 = L*S*R';

end