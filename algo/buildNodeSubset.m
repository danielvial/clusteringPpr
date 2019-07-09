% Sample nodes, either uniform or clustered (see Section 6.2)
%
% Input: G - properly formatted graph object (see README); k - number nodes
% to sample; type - sampling approach (uniform or clustered, denoted 'uni'
% or 'clulst'); clustType - clustering approach (1 or 2, differs in how
% nodes are sampled if building clustered set)tering type clustType (1 or 2)
%
% Output: U - sampled subset of nodes; condU - condutance of U, as defined
% in Section 6.2
%
function [U,condU] = buildNodeSubset(G,k,type,clustType)
    
    if isequal(type,'uni') % uniform sampling
        U = randsample(G.n,k);
    elseif isequal(type,'clust') % clustered sampling
        % for clustered sampling, we begin with a random node, then
        % greedily add nodes with probability proportional to number
        % (if clustType==1) or fraction (if clustType==2) of its incoming
        % edges originating in the subset; see Appendix H
        U = zeros(k,1); Sout = zeros(G.n,1); sampProb = Sout;
        U(1) = randi(G.n);
        for i=2:k
            % update sampling probabilities based on current subset
            Sout(G.Nout{U(i-1)}) = Sout(G.Nout{U(i-1)})+1;
            if clustType == 1
                sampProb(G.Nout{U(i-1)}) = Sout(G.Nout{U(i-1)});
            elseif clustType == 2
                sampProb(G.Nout{U(i-1)}) = Sout(G.Nout{U(i-1)})./G.Din(G.Nout{U(i-1)});
            end
            sampProb(U(1:i-1)) = 0; % don't sample node already in subset
            if nnz(sampProb) == 0
                % if sampling prob. zero for all nodes not currently in subset,
                % we've run out of nodes to add, so sample a uniform one
                newNode = randi(G.n);
                while ismember(newNode,U)
                    newNode = randi(G.n);
                end
            else
                % otherwise, add node with prob. prop. to sampling prob.
                newNode = randsample(G.n,1,true,sampProb);
            end
            U(i) = newNode;
        end
        % adding nodes in breadth-first-search fashion above could leave
        % unnatural artifacts, so we randomly permute to remove this effect
        U = U(randperm(length(U)));
    end
    
    % also compute conductance of the set (as described in paper,
    % conductance typically much lower when type='clust')
    condU = conductance(G,U);

end