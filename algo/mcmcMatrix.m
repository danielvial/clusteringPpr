% Joint multi-source walk sampling from Section 5.2 (matrix setting)
%
% Input: G - properly formatted graph object (see README); sigma - either
% the sigma_{avg} or sigma_{max} vector defined in Section 5.2; alpha -
% jump probability; w - number walks to sample
%
% Output: piHat - estimate of the matrix diag(1./sigma)*diag(sigma)*Pi
% defined in Section 5.2
%
function piHat = mcmcMatrix(G,sigma,alpha,w)

    w = ceil(w); % in case number walks not integer
    % sample starting nodes for walks; sigma is typically sparse in our
    % applications, faster to sample from only non-zero entries when using
    % MATLAB's randsample.m
    sigmaSupp = find(sigma);
    startV = sigmaSupp(randsample(length(sigmaSupp),w,true,full(sigma(sigmaSupp))));
    % sample ending locations for walks
    endV = walkSampler(G,startV,alpha);
    % construct output
    piHat = sparse(sigmaSupp,sigmaSupp,1./sigma(sigmaSupp),G.n,G.n)*...
        sparse(startV,endV,ones(w,1)/w,G.n,G.n);
    
end