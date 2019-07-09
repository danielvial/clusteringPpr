% Runs fwSingle.m for set of source nodes S
%
% Input: same as fwSingle.m, but here S is a set of source nodes
%
% Output: same as fwSingle.m, but here vectors stacked together as matrix
%
function [pS,rS,iter] = fwMultiSep(G,S,alpha,rmaxs)

    pS = sparse(G.n,length(S)); rS = pS; rS(S,:) = eye(length(S)); iter = 0;
    for i=1:length(S)
        [pS(:,i),rS(:,i),iterI] = fwSingle(G,S(i),alpha,rmaxs); iter = iter+iterI;
    end

end