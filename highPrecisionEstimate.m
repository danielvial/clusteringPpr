% Runs bwSingle.m to high precision to obtain (approximate) ground truth of
% true PPR values; used to validate parameters for multiple pair experiments
%
% Input: graphs - cell array of graph names; k - number of source/target
% pairs desired; alpha - jump probability; etaMult - bwSingle.m precision
% used here; deltaMult - precision to be used in singlePairComparison.m 
%
% Assumes data/graphs{i}.mat exists for each i (properly formatted graph
% object; see README)
%
% For each i, saves parameters and performance data to
% results/graphs{i}_hpe.mat
%
% See README for further explanation of how this function and
% singlePairComparison.m should be (jointly) used
%
function highPrecisionEstimate(graphs,k,alpha,etaMult,deltaMult)

    addpath algo; % sub-routines used
    
    for i=1:length(graphs) % loop over graphs
        
        load(['data/' graphs{i} '.mat']); % load graph data
        sigPairs = zeros(k,4); insigPairs = sigPairs; % initialize results
        
        for j=1:k % loop over source/target pairs
                        
            eta = etaMult/G.n; delta = deltaMult/G.n; % in paper, eta = 1/n, delta = 10/n
            t = randi(G.n); % choose random target
            
            % run bwSingle.m to high precision; sample one significant and
            % one insigificant source for this target (see Appendix A)
            [pt,rt] = bwSingle(G,t,alpha,eta);
            sig = find(pt>delta); sig = sig(randi(length(sig)));
            insig = find(pt<delta-eta); insig = insig(randi(length(insig)));
            
            % save source, target, pi_s(t) estimate, precision data
            sigPairs(j,:) = [sig,t,pt(sig),max(rt)]; insigPairs(j,:) = [insig,t,pt(insig),max(rt)];
            save(['results/' graphs{i} '_hpe.mat'],'sigPairs','insigPairs','etaMult','deltaMult');
            
        end
        
    end

end