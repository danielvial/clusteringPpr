% Estimates PPR vector via random walks
%
% Input: G - properly formatted graph object (see README); sigma - either a
% nonnegative vector whose length is number nodes in G (representing a
% distribution from which walks begin) or a scalar between 1 and number
% nodes in G (representing a source node from which walks begin); alpha -
% jump probability; wMult - number walks to samplenteger, take ceiling)
%
% Output: piHat - vector whose length is number nodes of v, estimate of
% pi_{sigma}(v), as defined in Section 2
%
function piHat = mcmcSingle(G,sigma,alpha,wMult)

    w = ceil(wMult); % in case number walks not integer
    
    % sample walks using walkSampler.m
    if length(sigma) == 1 % if given single node, all walks start from this node
        startV = sigma*ones(w,1);
        endV = walkSampler(G,startV,alpha);
    elseif length(sigma) == G.n % if given distribution, walks start at random states
        % sigma is typically sparse in our applications, faster to sample 
        % from only non-zero entries when using MATLAB's randsample.m
        sigmaSupp = find(sigma);
        startV = sigmaSupp(randsample(length(sigmaSupp),w,true,full(sigma(sigmaSupp))));
        endV = walkSampler(G,startV,alpha);
    else
        disp('error in mcmcSingle input'); piHat = []; return;
    end
    
    % construct output
    piHat = sparse(endV,ones(w,1),ones(w,1)/w,G.n,1);
    
end