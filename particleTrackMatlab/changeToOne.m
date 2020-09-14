function [pksUpdated,listUpdated] = changeToOne(npCombine,pksBest,list)
% objective: combines particles listed in 'npCombine' into a single
% centroid particle and removes these particle from 'list'

% list = complete list of dimers, trimers, etc...multimers
% pks = positions of particles in a single mulitmer

addpath('~/poincareProgs/particleTrackMatlab/');

pksUpdated = getCentroid(pksBest(npCombine,:));

npTemp = intersect(npCombine,list); %intersect(npCombine,list)) confirms that the particles listed in 'npCombine' are in 'list'
for jj = 1:numel(npTemp) %remove particles listed in npCombine from 'list' one at a time
    list(list == npTemp(jj)) = [];
end
listUpdated = list;
end