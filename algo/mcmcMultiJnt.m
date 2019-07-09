% Joint multi-source walk sampling from Section 5.1.1 (scalar setting)
%
% Input: G - properly formatted graph object (see README); rS - output from
% fwMultiSep.m; alpha - jump probability; wMult - vector whose length is
% number source nodes at hand (i.e. number columns of rS), wMult(i) =
% number random walk samples desired for i-th source for i-th source
%
% Output: piHat - matrix whose i-th column is estimate of pi_{sigma_s}, the
% PPR vector (see Section 2) for sigma_s = rs/||rs||_1 (see Section 4);
% totalW - number walks sampled (unlike mcmcMultiSep.m, not fixed a priori)
%
% Disclaimer: the implementation (specifically, the data structures used)
% is a bit odd and is not the obvious approach one would try after reading the
% paper; we found it's more efficient in practice than the obvious approach
%
function [piHat,totalW] = mcmcMultiJnt(G,rS,alpha,wMult)
    
    w = ceil(wMult); % in case number walks not integer
    cardS = size(rS,2); % total number source nodes
    
    % first, sample starting nodes for walks; data structures used here
    % are a bit odd, but essentially consist of three vectors/matrices:
    %   rowInd(i,j) = starting node for i-th walk of j-th source
    %   colInd(:,i) defined iteratively: colInd(1,i) = 1, and for j>1,
    %      if rowInd(j,i) = rowInd(j-1,i), colInd(j,i) = colInd(j-1,i)+1
    %      if rowInd(j,i) ~= rowInd(j-1,i), colInd(j,i) = 1
    %      (below, we don't actually use this iterative def'n, since we
    %      want to avoid loops, but iterative def'n more human-readable)
    %   maxCounts(i) will be max (across sources) number walks starting at i
    rowInd = zeros(max(w),cardS); colInd = rowInd; maxCounts = zeros(G.n,1); 
    for i=1:cardS
        % sample start nodes for walks for i-th source from corresponding
        % rs vector (see Section 4); rs is typically sparse, faster to
        % sample from only non-zero entries when using MATLAB's randsample.m
        rsSupp = find(rS(:,i));
        start = rsSupp(randsample(length(rsSupp),w(i),true,full(rS(rsSupp,i))));
        % next two lines basically create histogram of sampled start nodes;
        % histiNZ is list of nodes sampled at least once, histi is number
        % of times the corresponding nodes were sampled
        histi = sparse(start,ones(size(start)),ones(size(start)),G.n,1);
        histiNZ = find(histi); histi = histi(histiNZ);
        % update rowInd/colInd/maxCounts as described above (rowInd could've been
        % updated more easily from start vector used above, but we update from
        % histiNZ/hist instead since we need these vectors for the colInd update)
        rowInd(1:w(i),i) = repelem(histiNZ,histi);
        [colInd(1:w(i),i),~] = find(repmat(1:max(histi),[length(histi),1])'...
            <=repmat(histi,[1,max(histi)])');
        maxCounts(histiNZ) = max([maxCounts(histiNZ),histi],[],2);
    end
    
    % second, from each node, sample max number walks required by any source
    % (as in Section 5.1.1)
    totalW = sum(maxCounts); % total number walks we'll sample
    maxMaxCounts = max(maxCounts); % max number walks from any single node
    % next two lines construct startV, which lists each node from which we
    % sample walks the number of times we'll sample walks from it
    countsNZ = find(maxCounts);
    startV = repelem(countsNZ,maxCounts(countsNZ));
    [startVcnt,~] = find(repmat(1:maxMaxCounts,size(countsNZ))'...
        <=repmat(maxCounts(countsNZ),[1,maxMaxCounts])'); % will be used below
    endV = walkSampler(G,startV,alpha); % endV(i) = end node of walk from startV(i)
    
    % third, construct output as described above, using rowInd and colInd,
    % along with matrix walkMat, where walkMat(i,j) = end location of j-th
    % walk started from startV(i)
    walkMat = sparse(startV,startVcnt,endV,G.n,max(w));
    piHat = sparse(walkMat(sub2ind(size(walkMat),rowInd(rowInd>0),colInd(colInd>0))),...
        repelem(1:cardS,w),repelem(1./w,w),G.n,cardS);
    
end