% Compare FW-BW-MCMC and Bidirectional-PPR for single source/target pair;
% used to validate parameters for multiple pair experiments
%
% Inputs: graphs - cell array of graph names; alpha - jump probability; 
% allRmaxsO/allRmaxtO/allC0 - vector of rmaxs/rmaxt/c parameters for
% FW-BW-MCMC (each must have same size as graphs cell array);
% allRmaxtE/allCE - vector of rmaxt/c parameters for Bidirectional-PPR
% (each must have same size as graphs cell array)
%
% Assumes the existence of two files for each i:
%    data/graphs{i}.mat - properly formatted graph object (see README)
%    results/graphs{i}_hpe.mat - output of highPrecisionEstimate.m
% 
% For each i, saves parameters and performance data to
% results/graphs{i}_spc.mat
%
% See README for further explanation of how this function and
% highPrecisionEstimate.m should be (jointly) used
%
function singlePairComparison(graphs,alpha,allRmaxsO,allRmaxtO,allRmaxtE,allCO,allCE)

    addpath algo; % sub-routines used
    
    for i=1:length(graphs) % loop over graphs
        
        % load graph data
        load(['data/' graphs{i} '.mat']); load(['results/' graphs{i} '_hpe.mat']);
        
        % set parameters (note delta = 10/n for all experiments)
        delta = 10/G.n; cO = allCO(i); cE = allCE(i);
        rmaxsO = allRmaxsO(i); rmaxtO = allRmaxtO(i); rmaxtE = allRmaxtE(i);
        
        % initialize runtime/accuracy results
        times = zeros(4,2,length(sigPairs(:,1)),2); % time for each stage of algorithms
        error = zeros(2,length(sigPairs(:,1)),2); % error bounds (see Appendix A,H)
        
        for j=1:length(sigPairs(:,1)) % loop over source/target pairs

            % significant pairs, for which we want relative error bound (see Appendix A)
            s = sigPairs(j,1); t = sigPairs(j,2);
            % Bidirectional-PPR (only backward DP and MCMC; see Section 4)
            tic; [pt,rt] = bwSingle(G,t,alpha,rmaxtE); times(2,2,j,1) = toc;
            tic; piHat = mcmcSingle(G,s,alpha,cE*rmaxtE/delta); times(3,2,j,1) = toc;
            tic; exiEst = pt(s)+piHat'*rt; times(4,2,j,1) = toc;
            % FW-BW-MCMC (also uses forward DP; see Section 4)
            tic; [ps,rs] = fwSingle(G,s,alpha,rmaxsO); sumRs = sum(rs); times(1,1,j,1) = toc;
            tic; [pt,rt] = bwSingle(G,t,alpha,rmaxtO); times(2,1,j,1) = toc;
            tic; piHat = mcmcSingle(G,rs,alpha,cO*sumRs*rmaxtO/delta); times(3,1,j,1) = toc;
            tic; ourEst = pt(s)+ps'*rt+sumRs*piHat'*rt; times(4,1,j,1) = toc;
            % upper bound relative error (see Appendix A,H)
            error(1,j,1) = (abs(ourEst-sigPairs(j,3))+sigPairs(j,4))/sigPairs(j,3);
            error(2,j,1) = (abs(exiEst-sigPairs(j,3))+sigPairs(j,4))/sigPairs(j,3);
            
            % insignificant pairs, for which we want absolute error bound (see Appendix A)
            s = insigPairs(j,1); t = insigPairs(j,2);
            % Bidirectional-PPR (only backward DP and MCMC; see Section 4)
            tic; [pt,rt] = bwSingle(G,t,alpha,rmaxtE); times(2,2,j,2) = toc;
            tic; piHat = mcmcSingle(G,s,alpha,cE*rmaxtE/delta); times(3,2,j,2) = toc;
            tic; exiEst = pt(s)+piHat'*rt; times(4,2,j,2) = toc;
            % FW-BW-MCMC (also uses forward DP; see Section 4)
            tic; [ps,rs] = fwSingle(G,s,alpha,rmaxsO); sumRs = sum(rs); times(1,1,j,2) = toc;
            tic; [pt,rt] = bwSingle(G,t,alpha,rmaxtO); times(2,1,j,2) = toc;
            tic; piHat = mcmcSingle(G,rs,alpha,cO*sumRs*rmaxtO/delta); times(3,1,j,2) = toc;
            tic; ourEst = pt(s)+ps'*rt+sumRs*piHat'*rt; times(4,1,j,2) = toc;
            % check if error bound exceeds threshold (see Appendix A,H)
            error(1,j,2) = (abs(ourEst-insigPairs(j,3))+insigPairs(j,4))>2*exp(1)*delta;
            error(2,j,2) = (abs(exiEst-insigPairs(j,3))+insigPairs(j,4))>2*exp(1)*delta;
            
            % update results and save
            save(['results/' graphs{i} '_spc.mat'],...
                'delta','cO','cE','rmaxsO','rmaxtO','rmaxtE','times','error');
            
        end
        
    end

end