% Test performance of matrix approximation schemes from Section 5.2
% relative to baseline for both uniform and clustered source/target node
% subsets; used for Fig. 9 in Section 6.2
%
% Inputs: graphs - cell array of graph names; alpha - jump probability; 
% allRmaxsO/allRmaxtO/allC0 - vector of rmaxs/rmaxt/c parameters for
% proposed scheme (each must have same size as graphs cell array);
% allRmaxtE/allCE - vector of rmaxt/c parameters for baseline scheme
% (each must have same size as graphs cell array); trials - number of
% trials to run; l - number of sources/targets for each trial
%
% Assumes data/graphs{i}.mat exists for each i (properly formatted graph
% object; see README)
%
% For each i, saves parameters and performance data to
% results/graphs{i}_ma.mat
%
function matrixApprox(graphs,alpha,allRmaxsO,allRmaxtO,allRmaxtE,allCO,allCE,trials,l)

    addpath algo; % sub-routines used
    
    for i=1:length(graphs) % loop over graphs

        load(['data/' graphs{i} '.mat']); % load graph data

        % set parameters
        rmaxsO = allRmaxsO(i); rmaxtO = allRmaxtO(i); rmaxtE = allRmaxtE(i); 
        cO = allCO(i); cE = allCE(i);
        % for number walks, assume delta = 1/n
        wMultO = cO*rmaxtO/(10/G.n); wMultE = cE*rmaxtE/(10/G.n); 
        % to test accuracy, use high-precision estimate with tolerance 1/n        
        eta = 1/G.n; 

        % initialize performance results
        walks = zeros(2,3,trials); % number walks sampled
        times = walks; % runtime of algorithms
        clustQuant = zeros(2,2,trials); % clustering quantities from Section 5.2
        PiHatStBse = zeros(l,l,2,trials); % estimated matrix for baseline scheme
        PiHatStMax = PiHatStBse; % estimated matrix for sigma_{max} scheme (see Section 5.2)
        PiHatStAvg = PiHatStBse; % estimated matrix for sigma_{avg} scheme (see Section 5.2)
        PiStEta = PiHatStBse; % high-precision matrix estimate (for performance eval)
        etaPrac = zeros(2,trials); % eta is worse-case precision; record actual precision
        
        for j=1:trials % loop over trials

            for k=1:2 % k=1: uniform sources/targets, k=2: clustered sources/targets

                % sample sources, targets; also ensure surrogate matrix
                % discussed in Section 5.2 is nonzero
                surr = 0;
                while norm(surr,2) == 0
                    % sample sources and run forward DP
                    if k == 1
                        S = buildNodeSubset(G,l,'uni',0); % uniform node sampling
                    else
                        S = buildNodeSubset(G,l,'clust',1); % clustered node sampling
                    end
                    [pS,~] = fwMultiSep(G,S,alpha,rmaxsO);
                    % choose targets with large forward DP output (more
                    % likely that surrogate will be nonzero)
                    nonS = setdiff(1:G.n,S);
                    [~,idx] = sort(sum(pS(nonS,:),2),'descend');
                    T = nonS(idx(1:l));
                    % run backward DP and compute surrogate 
                    [pT,rT] = bwMultiSeq(G,T,alpha,rmaxtO);
                    surr = full(pT(S,:)+pS'*rT);
                end

                % baseline algorithm (backward DP + walks, see Section 5.2)
                tic;
                [pT,rT] = bwMultiSep(G,T,alpha,rmaxtE); % backward DP
                sigmaBse = sparse(S,1,1/l,G.n,1); % sigma matrix when no forward DP
                rSbse = sparse(S,1:l,1,G.n,l); % rS matrix when no forward DP
                wBse = wMultE*l; % number walks to sample
                PiHatStBse(:,:,k,j) = ... % sample walks and compute estimate
                    pT(S,:)+rSbse'*mcmcMatrix(G,sigmaBse,alpha,wBse)*rT;
                times(k,1,j) = toc;

                % proposed approaches (also uses forward DP, see Section 5.2)
                % forward/backward DP (used for both proposed approaches)
                tic;
                [pS,rS] = fwMultiSep(G,S,alpha,rmaxsO); % forward DP
                sumRs = sum(rS); % will use several times, so compute and save
                SigmaS = rS*sparse(1:l,1:l,1./sumRs,l,l); % sigma matrix
                [pT,rT] = bwMultiSeq(G,T,alpha,rmaxtO); % backward DP
                tmpTime = toc;
                times(k,2,j) = tmpTime; times(k,3,j) = tmpTime;
                % sigma_{max} approach (see Section 5.2)
                tic;
                sigmaMax = max(SigmaS,[],2); % un-normalized sigma_{max} vector
                l1inf = sum(sigmaMax); % clustering quantity/normalization factor
                sigmaMax = sigmaMax/l1inf; % normalized sigma_{max} vector
                wMax = wMultO*mean(sumRs)*l1inf; % number walks to sample
                PiHatStMax(:,:,k,j) = ... % sample walks and compute estimate
                    pT(S,:)+pS'*rT+rS'*mcmcMatrix(G,sigmaMax,alpha,wMax)*rT;
                times(k,2,j) = times(k,2,j)+toc;
                % sigma_{avg} approach (see Section 5.2)
                tic;
                surr = full(pT(S,:)+pS'*rT); % surrogate matrix
                srankSurr = (norm(surr,'fro')/norm(surr,2)); % clustering quantity
                sigmaAvg = mean(SigmaS,2); % sigma_{avg} vector
                wAvg = wMultO*mean(sumRs)*sqrt(l)*srankSurr; % number walks to sample
                PiHatStAvg(:,:,k,j) = ...% sample walks and compute estimate
                    pT(S,:)+pS'*rT+rS'*mcmcMatrix(G,sigmaAvg,alpha,wAvg)*rT;
                times(k,3,j) = times(k,3,j)+toc;

                % run backward DP to much higher precision (i.e. using eta
                % parameter) for accuracy evaluation
                [PiStEtaKJ,tmp] = bwMultiSep(G,T,alpha,eta); % backward DP
                PiStEta(:,:,k,j) = PiStEtaKJ(S,:); % sub-matrix of interest
                etaPrac(k,j) = max(sum(tmp,2))/l; % actual precision (see above)

                % update results and save
                walks(k,:,j) = [wBse,wMax,wAvg]; clustQuant(k,:,j) = [l1inf,srankSurr];
                save(['results/' graphs{i} '_ma.mat'],'S','T','walks','times','clustQuant',...
                    'PiHatStMax','PiHatStAvg','PiHatStBse','PiStEta','eta','l','etaPrac');

            end

        end

    end
    
end