% Approximate-Contributions algo by Andersen et al (Algorithm 2)
%
% Input: G - properly formatted graph object (see README); t - target node
% in G; alpha - jump probability; rmaxt - precision parameter
%
% Output: pt,rt - vectors defined in Section 4; iter - number iterations
% until algorithm converges (to evaluate performance)
%
function [pt,rt,iter] = bwSingle(G,t,alpha,rmaxt)

    pt = zeros(G.n,1); rt = zeros(G.n,1); rt(t) = 1; % initialization
    maxRt = 1; v = t; % initial choice of v^* in paper
    iter = 0; % number of iterations we've run
    
    % main loop for Algorithm 2
    while maxRt > rmaxt
        iter = iter+1;
        % update neighbors of v
        rt(G.Nin{v}) = rt(G.Nin{v})+(1-alpha)*rt(v)./G.Dout(G.Nin{v});
        % update v
        pt(v) = pt(v)+alpha*rt(v); rt(v) = 0;
        % choose next v
        [maxRt,v] = max(rt);
    end
    pt = sparse(pt); rt = sparse(rt); % return in sparse format since vectors now fixed
    
end