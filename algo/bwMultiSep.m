% Runs bwSingle.m (separately) for set of target nodes T
%
% Input: same as bwSingle.m, but here T is a set of target nodes
%
% Output: same as bwSingle.m, but here vectors stacked together as matrix
%
function [pT,rT,iter] = bwMultiSep(G,T,alpha,rmaxt)

    pT = sparse(G.n,length(T)); rT = pT; rT(T,:) = eye(length(T)); iter = 0;
    for i=1:length(T)
        [pT(:,i),rT(:,i),iterI] = bwSingle(G,T(i),alpha,rmaxt); iter = iter+iterI;
    end

end