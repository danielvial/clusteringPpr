% Creates directed, strongly-connected stochastic block model
%
% Input: n - number of nodes; c - number of communities; p -
% intra-community edge formation probability (each directed edge between
% nodes in same community independently present with probability p); q -
% inter-community edge formation probability (each directed edge between
% nodes in different communities independently present with probability q);
% suffix - optional string to append to output filename
%
% Typically will have p >> q to ensure densely-connected subsets of nodes
% (i.e. communities) that are sparsely inter-connected
%
% Assumes n is multiple of c (equal-sized communities, size n/c)
%
% Saves properly formatted graph (see README) to data/sbm_{suffix}.mat (if
% suffix provided) or data/sbm.mat (if no suffix provided)
%
function directConnectSbm(n,c,p,q,suffix)

    if nargin < 5
        savename = 'sbm.mat';
    else
        savename = ['sbm_' suffix '.mat'];
    end

    % will use David Gleich's gaimc (Graph Algorithms In Matlab Code)
    % package to ensure strong-connectedness; see README
    addpath gaimc;
    
    % construct SBM, continue until strongly-connected
    A = directSbm(n,c,p,q);
    % scomponents returns vector whose i-th entry is strongly-connected
    % component containing i-th node; when this vector contains one unique
    % value, the graph is strongly connected
    while length(unique(scomponents(A))) > 1
        A = directSbm(n,c,p,q);
    end
    
    % easier to construct adjacency matrix above, but now convert to format
    % used in experiments
    G.n = n; G.m = sum(A(:)); G.Dout = sum(A,2); G.Din = sum(A,1)';
    G.Nout = cell(n,1); G.Nin = cell(n,1);
    for i=1:n
        G.Nout{i} = find(A(i,:)); G.Nin{i} = find(A(:,i))';
    end
    
    % save graph
    save(savename,'G');
    
end

function A = directSbm(n,c,p,q)
    
    % first, sample intra-community edges (block diagonal matrix)
    A = binornd(1,p,n/c,n/c);
    for i=1:c-1
        A = blkdiag(A,binornd(1,p,n/c,n/c));
    end
    
    % next, sample inter-community edges (by creating all edges with
    % probability q each, then setting diagonal blocks to zero)
    B = binornd(1,q,n,n);
    for i=1:c
        B(1+(n/c)*(i-1):(n/c)*i,1+(n/c)*(i-1):(n/c)*i) = 0;
    end

    % finally, combine intra- and inter-community edges, and eliminate any
    % self-loops that may have formed
    A = A+B;
    A = A-diag(diag(A));

end