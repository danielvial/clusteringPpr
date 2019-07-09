% Cluster source nodes for distributed setting from Section 7
%
% Input: sigma_s - column-stochastic matrix (i-th col is normalized output
% of fwSingle.m or fwMulti.m for i-th source); k - number partitions
%
% Output: assign - vector whose length is number of sources (i.e. number
% cols of sigma_s) with assign(i) = partition i-th source assigned to
%
function assign = clusterSigma(sigma_s,k)
      
    [n,l] = size(sigma_s); % n = number nodes in graph, l = number sources
    assign = zeros(l,1); % initialize output described above
    
    % sigma_Si will contain centroids of clusters -- crucial that this is
    % sparse matrix! though MATLAB may complain the indexing below will be
    % slow, we found that this code runs hundreds of times faster when
    % using a sparse matrix versus a regular matrix (at least at the scale
    % of the experiments conducted for the paper)
    sigma_Si = sparse(n,k);
        
    % add random node to first cluster, update centroid
    init = randi(l); assign(init) = 1; sigma_Si(:,1) = sigma_s(:,init);
    
    % this for loop adds one node to each of remaining k-1 clusters,
    % randomly sampling in hopes of making the corresponding sigma vectors
    % far apart (similar to k-means++, see Section 7);
    % dist(i,j) will be distance between i-th source and j-th centroid
    dist = zeros(l,k-1);
    for i=2:k
        % update distances from sources to centroids
        dist(:,i-1) = sum(abs(sigma_s-repmat(sigma_Si(:,i-1),[1,l])))';
        % sample node for i-th partition proportional to distance
        minDist = min(dist(:,1:i-1),[],2); init = randsample(l,1,true,minDist);
        % add sampled node to i-th partition
        assign(init) = i; sigma_Si(:,i) = sigma_s(:,init);
    end
    % compute objective function value discussed in Section 7
    sigma_Si_l1 = full(sum(sigma_Si));
                    
    % iteratively add remaining nodes to partitions as in Section 7
    for i=1:l
        if assign(i) > 0
            continue; % do nothing if i already assigned
        end 
        % find partition jStar that, upon adding i, minimizes growth in
        % objective function value
        tmp = max(repmat(sigma_s(:,i),[1,k])-sigma_Si,0);
        d_s_Sj = sum(tmp); newL1inf = sigma_Si_l1+d_s_Sj;
        [L1inf_jStar,jStar] = min(newL1inf);
        % add i to partition jStar, update objective
        sigma_Si(:,jStar) = sigma_Si(:,jStar)+tmp(:,jStar);
        assign(i) = jStar; sigma_Si_l1(jStar) = L1inf_jStar;
    end
                
end


