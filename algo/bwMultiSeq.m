% Joint multi-target version of bwSingle.m from Section 5.1.2
%
% Input: same as bwSingle.m, but here T is a set of target nodes
%
% Output: same as bwSingle.m, but here vectors stacked together as matrix
%
function [pT,rT,iter,mc] = bwMultiSeq(G,T,alpha,rmaxt)

    pT = sparse(G.n,length(T)); rT = pT; rT(T,:) = eye(length(T)); % initialization
    iter = 0; % number of iterations we've run
    mc = 0; % number times we've used merge update from Section 5.1.2
    
    % iterate across target nodes T
    for i=1:length(T)
        pt = pT(:,i); rt = rT(:,i); maxRt = 1; v = T(i); % same as in bwMultiSep.m
        while maxRt > rmaxt
            iter = iter+1;
            tIdx = find(T(1:i-1)==v); % check if pv, rv known       
            if ~isempty(tIdx) % if pv, rv known, we can use merge update from Section 5.1.2
                mc = mc+1;
                rt(v) = 0; pt = pt+maxRt*pT(:,tIdx); rt = rt+maxRt*rT(:,tIdx);
            else % if pv, rv unknown, must use update from bwSingle.m
                rt(G.Nin{v}) = rt(G.Nin{v})+(1-alpha)*rt(v)./G.Dout(G.Nin{v});
                pt(v) = pt(v)+alpha*rt(v); rt(v) = 0;
            end
            [maxRt,v] = max(rt); % choose next v
        end
        pT(:,i) = pt; rT(:,i) = rt; % reduces number edits to sparse matrices (slow!)
    end

end