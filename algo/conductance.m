% Computes conductance of subset of nodes
%
% Input: Input: G - properly formatted graph object (see README); U -
% subset of nodes in G
%
% Output: cond - conductance of U, as defined in Section 6.1
%
function cond = conductance(G,U)

    % first compute numerator of expression in Section 6.1
    condNum = 0;
    for i=1:length(U)
        condNum = condNum+length(setdiff(G.Nout{U(i)},U));
    end
    
    % next compute denominator of expression in Section 6.1
    condDen = min([sum(G.Dout(U)),G.m-sum(G.Dout(U))]);
    
    % divide numerator by denominator to get condutance 
    cond = condNum/condDen;    
    
end