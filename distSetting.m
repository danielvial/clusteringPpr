% Compare heuristic and oracle methods against baseline for distributed
% setting experiments; used for Fig. 10 in Section 7
%
% Input: graphs - cell array of graph names; alpha - jump probability;
% allCO/allCE/allRmaxs/allRmaxt - vectors whose length are same length as
% graphs, allCO/allCE are c parameters for our scheme/baseline scheme
% (which use forward DP + walks/only walks), allRmaxs is forward DP
% precision parameter, allRmaxt is backward DP precision parameter (we
% assume backward DP was run offline to this tolerance); trials - number
% of trials to run; cardS - number of source nodes for each trial; clustNum
% - number partitions for clusterSigma.m; machNum - number machines
%
% Assumes cardS is multiple of clustNum, i.e. we can assign equal number
% sources to each machine
%
% Assumes data/graphs{i}.mat exists for each i (properly formatted graph
% object; see README)
%
% For each i, saves parameters and performance data to
% results/graphs{i}_ds.mat
%
function distSetting(graphs,alpha,allCO,allCE,allRmaxs,allRmaxt,trials,cardS,clustNum,machNum)

    addpath algo; % sub-routines used
    
    for i=1:length(graphs) % loop over graphs
        
        load(['data/' graphs{i} '.mat']); % load graph data
        
        % set parameters (note delta = 10/n for all experiments)
        delta = 10/G.n; cO = allCO(i); cE = allCE(i); 
        rmaxs = allRmaxs(i); rmaxt = allRmaxt(i);
        
        % initialize performance results
        times = zeros(3,3,trials); % time for each stage of algorithm and each method
        walks = zeros(3,trials); % number walks sampled for each method
        l1inf = zeros(2,trials); % objective function value for each method (see Section 7)

        for j=1:trials % loop over trials
            
            % sample sources; as in Section 7, we construct the source set
            % as disjoint union of clustNum clustered subsets; as shown in
            % paper, our heuristic method "recovers" these clusters
            S = zeros(1,cardS);
            for k=1:clustNum
                S(1+(cardS/clustNum)*(k-1):(cardS/clustNum)*k) ...
                    = buildNodeSubset(G,cardS/clustNum,'clust',1);
            end
     
            % baseline scheme described in Section 7 - arbitrarily assign
            % sources to machines (equal number to each machine),
            % independently sample walks for each source
            assign = ceil(randperm(cardS)/(cardS/machNum)); % arbitrary assignment
            % for simplicity, we don't actually run things in parallel, but
            % we simulate parallel performance by taking max time across
            % machines as runtime; for this, we use vector tmpTimes (also
            % tmpWalks, which lets us compute max number walks across machines
            tmpTimes = zeros(1,machNum); tmpWalks = tmpTimes;
            for k=1:machNum
                tic; 
                [~,tmpWalks(k)] = mcmcMultiSep(G,S(assign==k),alpha,cE*rmaxt/delta); 
                tmpTimes(k) = toc;
            end
            times(1,3,j) = max(tmpTimes); walks(1,j) = max(tmpWalks);
            
            % oracle scheme described in Section 7 - assign sources to
            % machines based on true clustering, run forward DP and sample
            % walks separately across machines
            assign = reshape(repmat(1:machNum,[(cardS/machNum),1]),[cardS,1]); % "true" assignment
            % in addition to tmpTimes, tmpWalks use above, here we'll also
            % compute objective function value discussed in Section 7
            tmpL1inf = zeros(1,machNum);
            for k=1:machNum
                % forward DP and walks for k-th machine
                tic;
                [~,rS] = fwMultiSep(G,S(assign==k),alpha,rmaxs);
                [~,tmpWalks(k)] = mcmcMultiJnt(G,rS,alpha,cO*sum(rS)*rmaxt/delta);
                tmpTimes(k) = toc;
                tmpL1inf(k) = sum(max(rS*sparse(1:size(rS,2),1:size(rS,2),1./sum(rS)),[],2));
            end
            times(2,3,j) = max(tmpTimes); walks(2,j) = max(tmpWalks); l1inf(1,j) = max(tmpL1inf);

            % heuristic scheme described in Section 7 - first assign
            % sources to machines arbitrarily and run forward DP; compute
            % heuristic partition based on forward DP; re-assign based on
            % this partition for walk sampling
            assign = ceil(randperm(cardS)/(cardS/machNum)); % arbitrary assign for forward DP
            rS = sparse(G.n,cardS); % will need rs vectors for every source to do walk partition
            for k=1:machNum
                % forward DP for k-th machine
                tic;
                [~,rS(:,assign==k)] = fwMultiSep(G,S(assign==k),alpha,rmaxs);
                tmpTimes(k) = toc;
            end
            times(3,1,j) = max(tmpTimes);
            % compute assignment using heuristic partitioning scheme
            tic; assign = clusterSigma(rS*sparse(1:size(rS,2),1:size(rS,2),1./sum(rS)),machNum); 
            times(3,2,j) = toc;
            for k=1:machNum
                % walks for k-th machine
                tic;
                [~,tmpWalks(k)] = mcmcMultiJnt(G,rS(:,assign==k),alpha,...
                    cO*sum(rS(:,assign==k))*rmaxt/delta);
                tmpTimes(k) = toc;
                rSk = rS(:,assign==k);
                tmpL1inf(k) = sum(max(rSk*sparse(1:size(rSk,2),1:size(rSk,2),1./sum(rSk)),[],2));
            end
            times(3,3,j) = max(tmpTimes); walks(3,j) = max(tmpWalks); l1inf(2,j) = max(tmpL1inf);

            % update results and save
            save(['results/' graphs{i} '_ds.mat'],'times','walks','l1inf',...
                'delta','cO','cE','rmaxs','rmaxt','trials','cardS','clustNum','machNum');
        
        end
        
    end

end



