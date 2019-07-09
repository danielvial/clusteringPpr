% Creates directed, strongly-connected Erdos Renyi graph
%
% Input: n - number of nodes; p - edge formation probability (each directed
% edge will be present with probability p, independent across edges);
% suffix - optional string to append to output filename
%
% Saves properly formatted graph (see README) to data/er_{suffix}.mat (if
% suffix provided) or data/er.mat (if no suffix provided)
%
function directConnectEr(n,p,suffix)
    
    if nargin < 3
        savename = 'er.mat';
    else
        savename = ['er_' suffix '.mat'];
    end

    % will use David Gleich's gaimc (Graph Algorithms In Matlab Code)
    % package to ensure strong-connectedness; see README
    addpath gaimc; 
    
    % construct ER graph, continue until strongly-connected
    A = directEr(n,p);
    % scomponents returns vector whose i-th entry is strongly-connected
    % component containing i-th node; when this vector contains one unique
    % value, the graph is strongly connected
    while length(unique(scomponents(A))) > 1
        A = directEr(n,p);
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

function A = directEr(n,p)

    A = binornd(1,p,n,n); % each directed edge independently present with prob. p
    A = A-diag(diag(A)); % eliminate any self-loops that formed
    
end