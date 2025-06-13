function [L,R] = nmf2qdr(N)
%NMF2QDR Summary of this function goes here
%   Detailed explanation goes here

[L0,R0] = qdrinit(N);

[L,R] = ANLS(N,L0,R0);


end

