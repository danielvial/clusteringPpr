% Compare FW-BW-MCMC and Bidirectional-PPR for multiple source/target pairs
% (scalar estimation setting); used for Fig. 7-8 in Section 6.2
%
% Inputs: graphs - cell array of graph names; alpha - jump probability; 
% allRmaxsO/allRmaxtO/allC0 - vector of rmaxs/rmaxt/c parameters for
% FW-BW-MCMC (each must have same size as graphs cell array);
% allRmaxtE/allCE - vector of rmaxt/c parameters for Bidirectional-PPR
% (each must have same size as graphs cell array); trials - number of
% trials to run; k - number of sources/targets for each trial; type - 'uni'
% or 'clust', clustType - 1 or 2 (see algo/buildNodeSubset)
%
% Assumes data/graphs{i}.mat exists for each i (properly formatted graph
% object; see README)
%
% For each i, saves parameters and performance data to
% results/graphs{i}_mpc.mat
%
function multiPairComparison(graphs,alpha,allRmaxsO,allRmaxtO,allRmaxtE,allCO,allCE,...
    trials,k,type,clustType)

    addpath algo; % sub-routines used
    
    for i=1:length(graphs) % loop over graphs
        
        load(['data/' graphs{i} '.mat']); % load graph data
        
        % set parameters (note delta = 10/n for all experiments)
        delta = 10/G.n; cO = allCO(i); cE = allCE(i);
        rmaxsO = allRmaxsO(i); rmaxtO = allRmaxtO(i); rmaxtE = allRmaxtE(i);
                
        % initialize performance results
        cond = zeros(2,trials); % conductance of sources/targets chosen
        times = zeros(4,2,trials); % time for each stage of algorithms
        % clustering quantities (see Section 5)
        clustS = zeros(1,trials); % source quantity
        clustT = zeros(1,trials); % target quantity
        clustST = zeros(2,trials); % joint source/target quantity (stable rank)
        % algorithmic performance details (see Sections 4-5)
        bwIter = zeros(2,trials); % number backward DP iterations
        mc = zeros(1,trials); % number times merge update used
        fwIter = zeros(1,trials); % number forward DP iterations
        walks = zeros(2,trials); % number random walks sampled
        
        for j=1:trials % loop over trials
            
            % sample sources/targets
            [S,cond(1,j)] = buildNodeSubset(G,k,type,clustType);
            [T,cond(2,j)] = buildNodeSubset(G,k,type,clustType);
            
            % existing method (see start of Section 5.1)
            tic; [pT,rT,bwIter(2,j)] = bwMultiSep(G,T,alpha,rmaxtE); times(2,2,j) = toc;
            tic; [piSigHat,walks(2,j)] = mcmcMultiSep(G,S,alpha,cE*rmaxtE/delta); times(3,2,j) = toc;
            tic; piHatStE = pT(S,:)+piSigHat'*rT; times(4,2,j) = toc;
            
            % our method (see Section 5)
            tic; [pS,rS,fwIter(j)] = fwMultiSep(G,S,alpha,rmaxsO); sumRs = sum(rS); times(1,1,j) = toc;
            tic; [pT,rT,bwIter(1,j),mc(j)] = bwMultiSeq(G,T,alpha,rmaxtO); times(2,1,j) = toc;
            tic; [piSigHat,walks(1,j)] = mcmcMultiJnt(G,rS,alpha,cO*sumRs*rmaxtO/delta); times(3,1,j) = toc;
            tic; piHatStO = pT(S,:)+pS'*rT+diag(sumRs)*piSigHat'*rT; times(4,1,j) = toc;
            
            % compute clustering quantities (since actual PPR matrix unknown,
            % can only estimate stable rank; do so using both estimates
            clustS(j) = sum(max(rS*diag(1./sum(rS)),[],2));
            clustT(j) = sum(sum(tril(pT(T,:),-1)>rmaxtO));
            clustST(1,j) = (norm(piHatStO,'fro')/norm(full(piHatStO),2))^2;
            clustST(2,j) = (norm(piHatStE,'fro')/norm(full(piHatStE),2))^2;
            
            % update results and save
            save(['results/' graphs{i} '_' type '_' num2str(clustType) '_mpc.mat'],...
                'delta','cO','cE','rmaxsO','rmaxtO','rmaxtE',...
                'cond','bwIter','mc','walks','fwIter',...
                'times','clustS','clustT','clustST');
            
        end
        
    end

end



