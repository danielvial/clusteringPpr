% Samples random walks for use in PPR estimation
%
% Input: G - properly formatted graph object (see README); startV -
% locations from which walks begin; alpha - jump probability
%
% Output: endV - vector with same dimension as startV, endV(i) is node at
% which i-th walk ends
%
function endV = walkSampler(G,startV,alpha)
    
    w = length(startV); % number walks
    L = geornd(alpha,[w,1]); % walk lengths
    endV = zeros(w,1); % end locations of walks
    
    % iteratively sample walks
    for i=1:w
        loc = startV(i);
        for j=1:L(i) % sample next state from neighbors of current state
            NoutLoc = G.Nout{loc}; loc = NoutLoc(randi(length(NoutLoc)));
        end
        endV(i) = loc;
    end

end