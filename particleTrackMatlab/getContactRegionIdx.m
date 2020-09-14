function [contactRegionIdx,boolMaxima] = getContactRegionIdx(profileGreen,centerIdx,buffer)
addpath('~/poincareProgs/particleTrackMatlab/');

%objective: to return the region between particles' edges
% MAJOR Assumption: the local maximium closest to centerIdx is the locates
%                   a particle's edge

% contactRegionIdx = indices of the local maxima in profileGreen that
%                    identify the edges of particles
% 
% boolMaxima(1) = 1 if contactRegionIdx(1) is a TRUE local maxima returned
%                 by getLocalMaxima; = 0 if contactRegionIdx(1) is just an
%                 extraploted location
% similarly, boolMaxima(2) = 1 if contactRegionIdx(2) is a TRUE local
% max. otherwise boolMaxima(2) == 0

% profileGreen = linear intensity profile between 2 pts in InormGreen
% centerIdx = the center profileGreen

% Example:
% profileGreen is:
%                           contactPt
%                              or
%    pt1           (+)<-----centerIdx ------>(-)                     pt2
%                              |
%                              |
% (start)--x-------x--------(center)---------------x---x----------(stop)
%      1   2   3   4   5   6   7   8   9   10  11  12  13   14  15  17
% 
% localMaximaIdx = [2,4,12,13]; %max peaks returned by getLocalMaxima(profileGreen)
% centerIdx = 7
% 
% displacement = localMaximaIdx - centerIdx 
%             = [2,4,12,13] - 7
%             = [-5,-3,5,6]
% 
% those peaks closest to centerIdx provide the most conservative estimates
% for the edges of the particle. 
% 
% nanmin(abs(displacement(displacement<=0))) %smallest displacement to the LEFT of centerIdx
% nanmin(abs(displacement(displacement>=0))) %smallest displacement to the RIGHT of centerIdx
% 
% contactRegionStart = centerIdx - nanmin(abs(displacement(displacement<=0))); %index to the left of centerIdx
% contactRegionStop = centerIdx + nanmin(abs(displacement(displacement>=0))); %index to the left of centerIdx

if(~exist('buffer','var'))
    buffer = 0;
end

boolMaxima = zeros([1,2]);

localMaximaIdx = getLocalMaximaB(profileGreen);

if(~isempty(localMaximaIdx))
    displacement = localMaximaIdx - centerIdx;

    if(numel(displacement(displacement<=0))>0) %index to the LEFT of centerIdx
        contactRegionStart = centerIdx - nanmin(abs(displacement(displacement<=0)));
        boolMaxima(1) = 1;
    end
        
    if(numel(displacement(displacement>=0))>0) %index to the RIGHT of centerIdx
        contactRegionStop = centerIdx + nanmin(abs(displacement(displacement>=0)));
        boolMaxima(2) = 1;
    end
    
    if(isempty(displacement(displacement<=0))) %There's no localMaximaIdx to the LEFT of centerIdx
        contactRegionStart = centerIdx - nanmin(abs(displacement(displacement>=0))); %location on profileGreen that is symmetric wrt contactRegionStop BUT is NOT a local maxima
    end
    
    if(isempty(displacement(displacement>=0))) %There's no localMaximaIdx to the RIGHT of centerIdx
        contactRegionStop = centerIdx + nanmin(abs(displacement(displacement<=0)));
    end
    
    %add buffer...
    contactRegionStart = contactRegionStart - buffer;
    contactRegionStop = contactRegionStop + buffer;
    
    %correct in case 'buffer' changes contactRegionStart/Stop to unrealistic values
    if(contactRegionStart < 1)
        contactRegionStart = 1;
    end
    
    if(contactRegionStop > numel(profileGreen))
        contactRegionStop = numel(profileGreen);
    end
end

if(isempty(localMaximaIdx)) %in the event getLocalMaxima() is unable to find a max peak in profileGreen
    contactRegionStart = 1;
    contactRegionStop = numel(profileGreen);
end

contactRegionIdx = [contactRegionStart,contactRegionStop];
end