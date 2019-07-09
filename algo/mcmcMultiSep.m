% Runs mcmcSingle.m (separately) for set of source nodes S
%
% Input: same as mcmcSingle.m, but here S is a set of source nodes
%
% Output: same as mcmcSingle.m, but here  vectors stacked together as
% matrix; also report total number of walks sampled (this may seem silly,
% since number walks is also an input, but this way output has same format
% as mcmcMultiJnt.m, where number walks not fixed a priori)
%
function [piHat,totalW] = mcmcMultiSep(G,S,alpha,wMult)

    w = ceil(wMult); % in case number walks not integer
    cardS = length(S); totalW = w*cardS; % total number walks we'll sample
    
    % format input for walkSampler.m
    startV = repmat(reshape(S,[1,length(S)]),[1,w]); 
    startVind = repmat(1:cardS,[1,w]);
    endV = walkSampler(G,startV,alpha);
    
    % construct output
    piHat = sparse(endV,startVind,ones(size(startV))/w,G.n,cardS);

end