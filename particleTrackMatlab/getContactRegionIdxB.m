function [contactRegionIdx,contactIdx] = getContactRegionIdxB(profileGreen)
addpath('~/poincareProgs/particleTrackMatlab/');

% Author: Eru K.
% Date: 31-May-2016
% adopted from getContactRegionIdxB.m


% StartContactStopIdx =
% [contactRegionIdx(1),contactIdx,contactRegionIdx(2)]
% 
% contactRegionIdx(1) = largest peak in profileGreen locating edge of
%                       particle np
% 
% contactRegionIdx(2) = largest peak in profileGreen locating edge of
%                       particle npp, which neighbors particle np
% 
% contactIdx == location of local minimum between the 2 largest peaks
%               in profileGreen. The is EXPERIMENTALLY and DEFINITIVELY
%               locates the space between the neighboring particles

localMaximaIdx = getLocalMaximaB(profileGreen);


%there should ALWAYS be at least one peak in profileGreen that locating the
%edges of neighboring particles.  In the event that there are no peaks,
%return the entire profileGreen as the contactRegion.  contactIdx == nan
if(isempty(localMaximaIdx))
    contactRegionIdx(1) = 1;
    contactRegionIdx(2) = numel(profileGreen);
    contactIdx = 1; %CanNOT let this be 'nan' because MATLAB crashes when attempt to locate and index nan
end

if(numel(localMaximaIdx)==1)%Only one localMaxima. Replace defualt contactRegionIdx(1) and contactRegionIdx(2) with the single localMaxima
    contactRegionIdx(1) = localMaximaIdx;
    contactRegionIdx(2) = localMaximaIdx;
    contactIdx = localMaximaIdx;
end

if(numel(localMaximaIdx)>1)
    [naught,idxTemp] = sort(profileGreen(localMaximaIdx),'descend'); %arrange localMaximaIdx from largest to smallest values of profileGreen(localMaximaIdx)
    localMaximaIdxSorted = localMaximaIdx(idxTemp);

    contactRegionIdx(1) = nanmin([localMaximaIdxSorted(1),localMaximaIdxSorted(2)]);%edge of particle np
    contactRegionIdx(2) = nanmax([localMaximaIdxSorted(1),localMaximaIdxSorted(2)]);%edge of neighboring particle is farther away (i.e., larger index on profileGreen)

    %locate the local minimum.  This is the middle of the contact region!
    contactTempIdx = contactRegionIdx(1):contactRegionIdx(2);
    contactIdx = nanmin(contactTempIdx(profileGreen(contactTempIdx) == nanmin(profileGreen(contactTempIdx)))); %local minimum between contactRegionIdx(1) and contactRegionIdx(2). nanmin(...) returns ONLY one index in rare event that there 2 indices in contactTempIdx where profileGreen(contactTempIdx) == nanmin(profileGreen(contactTempIdx))
end
end