function [idxMin,idxContact,idxMax,contactWidthIdx,boolContact] = getContactContrastIdx(posR,profileGreen,profileRed,contactRegionIdx)
% Author: Eru K.
% Date: 26-May-2016
% objective: return indices inside contactRegion 
% 
%   when boolContact == 1, it returns all indices inside the contactRegion
%   where profileGreen > profileRed 
% 
%   when boolContact == 0, it returns all indices inside the contactRegion
%   where profileGreen <= profileRed 
% 
% boolContact == 1 when there's at least one index inside the contactRegion
% where profileGreen > profileRed

tempDiff = profileGreen-profileRed;
contrastGreenRed = abs(profileGreen-profileRed);
idxContrast = 1:numel(contrastGreenRed);

if(numel(intersect(find(tempDiff>0),contactRegionIdx(1):contactRegionIdx(2)))>0)%there's at least one locating inside the contact region where profileGreen > profileRed
    boolContact = 1;
    contactWidthIdx = intersect(find(tempDiff>0),contactRegionIdx(1):contactRegionIdx(2));

    idxContrast(contactWidthIdx) = [];
    contrastGreenRed(idxContrast) = [];%deletes all indices except contactWidthIdx
    
    %contactWidthIdx and contrastGreenRed are the same length and the represent the same locations
    idxMin = nanmin(contactWidthIdx(contrastGreenRed == nanmin(contrastGreenRed)));%location of the smallest value in contrastGreenRed
    idxMax = nanmin(contactWidthIdx(contrastGreenRed == nanmax(contrastGreenRed)));

    idxWrtContact = abs(contactWidthIdx - find(posR==0));
    idxContact = contactWidthIdx(idxWrtContact == nanmin(idxWrtContact)); %numel(idxContact) == 1 or 2, e.g., [1 3 4 4 5] ==> idxContact = 2,3  
end

if(numel(intersect(find(tempDiff>0),contactRegionIdx(1):contactRegionIdx(2)))<=0)%there's no value where profileGreen > profileRed
    boolContact = 0;
    contactWidthIdx = intersect(find(tempDiff<=0),contactRegionIdx(1):contactRegionIdx(2));
    
    idxContrast(contactWidthIdx) = [];
    contrastGreenRed(idxContrast) = [];%deletes all indices except contactWidthIdx

    %contactWidthIdx and contrastGreenRed are the same length and the represent the same locations
    idxMin = nanmin(contactWidthIdx(contrastGreenRed == nanmin(contrastGreenRed)));%location of the smallest value in contrastGreenRed
    idxMax = nanmin(contactWidthIdx(contrastGreenRed == nanmax(contrastGreenRed)));

    idxWrtContact = abs(contactWidthIdx - find(posR==0));
    idxContact = contactWidthIdx(idxWrtContact == nanmin(idxWrtContact));
end
end

