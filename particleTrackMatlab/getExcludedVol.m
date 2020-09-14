function [excludedVol,rSmall,rLarge,dGap] = getExcludedVol(profile,edgesIdx)
% Author: Eru K.
% objective: Cacluate excluded volume

% input variables self explanatory:
%   profile = linear intensity profile between 2 particles (call them
%             particles 1 to 2)
% 
%   contactRegion = indices along profile that identifies inter-particle
%                   space between particles 1 and 2
% 
%   edgesIdx = 1 or 2 largest of the local maxima in 'profile' idenifying
%   particles' edges
%       numel(egdes) == 1 when particles are close 
%       numel(egdes) == 2 when particles are far apart 

volExclude = nan([3,1]);
dGap = nan(1);

if(numel(edgesIdx) == 1)%only 1 peak in 'profile' so particles are very close maybe touching
    rSmall = min(numel(profile)-edgesIdx,edgesIdx-1); %minus 1 because you need distance and MATLAB starts indexing at 1 instead of 0
    rLarge = max(numel(profile)-edgesIdx,edgesIdx-1);
end

if(numel(edgesIdx) == 2) %2 peaks in 'profile'. particles' edges are far apart.  Therefore, find the edges
    distFromStart = min(edgesIdx(1)-1,edgesIdx(2)-1);
    distFromEnd = min(numel(profile)-edgesIdx(1),numel(profile)-edgesIdx(2));

    rSmall = min(distFromStart,distFromEnd); %define smallest distance as radius of sphere1
    rLarge = max(distFromStart,distFromEnd); %assign larger distance to radius of sphere2
    dGap = numel(profile) - rSmall - rLarge; %inter-particle space
end


% |       |     |     |\
% |rSmall |     |     | \rLarge
% |       |     |    a|  \
% |       |     |     |   \
%  -------|     |-----|----x
%         |     |  h  | b /
% sphere1 |inter-    sphere2
%         |space|     | /
%         |     |     |/

volExclude(1) = (rSmall^3)*(4 - (2*pi/3));%half the volume of the cube surrounding the half of spherical partice1 

%particle2 is the same size or larger than particle1. If larger, then
%only a spherical cap of particle 2 is in the exlcuded volume
% 
% rLarge = radius of sphere2
% b = distance for center of sphere2 and ciruclar base spherical cap
% a = radius of spherical cap's circular base (a = radius of sphere1)
% h = height of spherical cap
%  
%      |\
%      | \rLarge
%     a|  \
%      |   \
%  ----|----x
%   h  | b /
%      |  /
%      | /
%      |/
%
% Pythagoras: rLarge^2 = a^2 + b^2
%             b = sqrt(rLarge^2 - a^2)
% 
% h + b = rLarge
%     h = rLarge - b 
%     h = rLarge - sqrt(rLarge^2 - a^2)
%     h = rLarge - sqrt(rLarge^2 - rSmall^2)
h = rLarge - sqrt((rLarge^2) - (rSmall^2));
volExclude(2) = (4*h*(rSmall^2)) - (pi*h/6)*((3*(rSmall^2))+(h^2));

volExclude(3) = 4*dGap*(rSmall^2);

excludedVol = nansum(volExclude(:));
end
