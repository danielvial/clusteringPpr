% Approximate-PageRank algo by Andersen et al (Algorithm 1)
%
% Input: G - properly formatted graph object (see README); s - source node
% in G; alpha - jump probability; rmaxs - precision parameter
%
% Output: ps,rs - vectors defined in Section 4; iter - number iterations
% until algorithm converges (to evaluate performance)
%
function [ps,rs,iter] = fwSingle(G,s,alpha,rmaxs)

    ps = zeros(G.n,1); rs = zeros(G.n,1); rs(s) = 1; % initialization
    dInvRs = rs/G.Dout(s); % dedicated rs./Dout vector for v^* as in Algorithm 1
    maxDinvRs = dInvRs(s); v = s; % initial choice of v^*
    iter = 0; % number of iterations we've run
    
    % main loop for Algorithm 1
    while maxDinvRs > rmaxs
        iter = iter+1;
        % update neighbors of v
        rs(G.Nout{v}) = rs(G.Nout{v})+(1-alpha)*rs(v)/G.Dout(v);
        dInvRs(G.Nout{v}) = rs(G.Nout{v})./G.Dout(G.Nout{v});
        % update v
        ps(v) = ps(v)+alpha*rs(v); rs(v) = 0; dInvRs(v) = 0;
        % choose next v
        [maxDinvRs,v] = max(dInvRs);
    end
    
    ps = sparse(ps); rs = sparse(rs); % return in sparse format since vectors now fixed
    
end