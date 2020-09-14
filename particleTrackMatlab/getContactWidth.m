function [boolContactPt,contactWidthIdx] = getContactWidth(posR,contactRegionIdx,profileGreen,profileRed)
% Author: Eru K.
% Date: 11-June-2016
% objective: to return indices defining the 'breadth' of the contactRegion
% Important note: the breadth of the contactRegion must satisfy BOTH
% criteria of fluorphore exclusion
%   1) one peak in profileGreen(posR(ii)). 
%       i.e., contactRegion(posR(ii)) = 0
%   2) profileGreen(posR(ii)) > profileRed(posR(ii))


%                           ..........................contactWidthIdx(1)
%                         |
%                        (+)   ...
%                         |  posR(contactIdx-2) = -2
%                         |  posR(contactIdx-1) = -1
%  pksBest(np,:)---------posR(contactIdx) = 0---------------pksBest(npp,:)
%                         |  posR(contactIdx+1) = 1
%                         |  posR(contactIdx+2) = 2
%                        (-)   ...
%                         |............................contactWidthIdx(2)
%                   |<----|---->| (contactRegion = 8)
%                  |<-----|----->| (contactRegion = 10)
%                 |<------|------>| (contactRegion = 12)

% contactRegion = distacne btwn peaks in profileGreen(*) along the line
%                 btwn pksBest(np,:) and pksBest(npp,:)
% 
% (*) == this profileGreen is different from 'profileGreen' measured along
%        the perpendicular line btwn pksBest(np,:) and pksBest(npp,:)

% contactRegionIdx = [startIdx, contactIdxGeo, stopIdx]

% --------------------------x
% Step1: preliminary check to see that profileGreen/Red are measured
% between particles pksBest(np,:) and pksBest(npp,:) to begin with
boolInside = nan(size(posR));
% boolInside == nan for single peaks in profileGreen, 
% boolInside == 1 if there are multiple peaks in profileGreen and
%               contactIdxXY is located btwn the 2 largest peaks, which are
%               presumably the particles' edges  
% boolInside == 0 if there's multiple peaks in profileGreen and
%               contactIdxXY is NOT located between the 2 largest peaks
for ii = 1:numel(posR)  
    if((contactRegionIdx(ii,3) - contactRegionIdx(ii,1)) > 0) %non-zero/non-vanishing contactRegion. This means there are 2 peaks in profileGreen
        if((contactRegionIdx(ii,2) >= contactRegionIdx(ii,1)) && (contactRegionIdx(ii,2) <= contactRegionIdx(ii,3)))%verifies contactIdx is located inside the contactRegion btwn contactRegionIdx(1) and contactRegionIdx(2)
            boolInside(ii,1) = 1;%non-vanishing contactRegion where the contact pt. [i.e., contactRegionIdx(ii,2)] is located somewhere INSIDE the contactRegion
        else
            boolInside(ii,1) = 0;%non-vanishing contactRegion where the contact pt. [i.e., contactRegionIdx(ii,2)] is located somewhere OUTSIDE the contactRegion
        end
    end
end


% Size of the contactRegion.
% What's the distance btwn contactRegionIdx(ii,3) and
% contactRegionIdx(ii,1)
% 
% Validity of the contactRegion.
% IS contactRegionIdx(ii,2) btwn 2 peaks in profileGreen [i.e.,
% contactRegionIdx(ii,1) and contactRegionIdx(ii,3)]?
% 
% Fluorescent signals inside the contact region.
% Finally, what's the value of profileRed and profileGreen at different
% locations between contactRegionIdx(ii,1) and contactRegionIdx(ii,3)

% ----------------------x
% Things to do:
% Need to see how boolStart/Stop are determined in:
% 04d_distWrtContact_contactWidthStartIdx_contactIdx_StopIdx_boolStart_boolStop_GREENmax_contact_min_REDmax_contact_min
% ...because you're going to use these to determine validity of the
% contactRegion. 
%   contactRegionSize(ii) ~= contactRegionIdx(ii,3) - contactRegionIdx(ii,1) 
%   if boolStart == 0 or boolStop == 0
% 
% ANSWER: boolStart/Stop are important!
% if boolStart == 1, then contactRegionIdx(ii,1) is an actual peak in
%       profileGreen and therefore an actual particle's edge
% 
% otherwise, boolStart == 0 means contactRegionIdx(ii,1) is NOT a
% particle's fluorescent edge
% 
% This is important because 
%   contactRegion = contactRegionIdx(ii,3) - contactRegionIdx(ii,1) 
%   if either contactRegionIdx(ii,3) or contactRegionIdx(ii,1) is NOT real,
%   then you can NOT accurately calculate the size of the contact region.


% NEXT,
%   plot [profileGreen/Red vs. posR]
%   AND [posR vs. contactRegion] for profileGreen/Red in XY, XZ, and YZ
%   directions

% ------------------x

% Step1b: if most contactRegionIdx(ii,2) has been confirmed to be located
% inside the contact region for most values of ii...
% 
%   i.e., contactRegionIdx(ii,2) >= contactRegionIdx(ii,1) 
%   and 
%   contactRegionIdx(ii,2) <= contactRegionIdx(ii,3)
% 
% ...then locate those indices where profileGreen(ii) >  profileRed(ii)
idxTouching = 1:numel(posR);
contrastGreenRedtemp = profileGreen - profileRed;
if(numel(boolInside(boolInside==1)) > numel(boolInside(boolInside==0)))%verifies contactRegionIdx(ii,2) is inside the contactRegion for most values of ii
    idxTouch = intersect(idxTouching(contrastGreenRedtemp > 0),idxTouching(isnan(boolInside) == 1)); %indices (aka values of ii) where profileGreen(ii) > profileRed(ii) AND where there's a single peak in profileGreen(ii) [i.e., boolInside(ii) == nan]
    diffTemp = idxTouch - idxTouching(posR==0);
    % -----(-)-----<-----------------|--------------->---------(+)-----
    %   posR < 0                  posR = 0                  posR > 0
else
    diffTemp = [];
end

% step2: return boolContactPt and contactWidthIdx 
% isnan(diffTemp)==1 means there are not any values of ii
% where:  
%   1) profileGreen(ii) > profileRed(ii) 
%   AND
%   2) there's a single peak in profileGreen
% 
% % % % % % 
% experimental evidence that
%                    profileGreen/Red measures fluorescence btwn
%                    pksBest(np,:) and pksBest(npp,:) [and orthogonal to
%                    the contactPt]   
% 
% (isempty(diffTemp)==1) means intersect(...,...) could NOT find any single
%                        peaks in profileGreen(posR(ii)) that are greater
%                        than profileRed(posR(ii))
if(isempty(diffTemp)==1)
    boolContactPt = 0;%NO, contactPt is NOT located around tempXY/XZ/YZ(1,:) == 0
    contactWidthIdx = [nan,nan];
end

if((numel(diffTemp(diffTemp<=0)) > 0)  && (numel(diffTemp(diffTemp>=0)) > 0))%checks to see if there are indices to the right AND left of tempXY/XZ/YZ(1,:)==0 with a single peak in profileGreen that is > profileRed
    idxDiff = 1:numel(diffTemp);
    contactWidthIdxStart = idxTouch(diffTemp == nanmax(diffTemp(diffTemp<=0)));
    contactWidthIdxStop = idxTouch(diffTemp == nanmin(diffTemp(diffTemp>=0)));
    if(contactWidthIdxStart == contactWidthIdxStop)%ONLY occurs when contactWidthIdxStart == 0 and contactWidthIdxStop == 0

        %        contactWidthIdxStart   
        %                 |
        % -----(-)<-------0------>(+)-----
        %                 |     
        %         contactWidthIdxStop

        boolContactPt = 1;%YES, contactPt is located around tempXY/XZ/YZ(1,:) == 0

        % 1st condition makes sure kk is NEVER > numel(idxTouch); 
        % 2nd condition requires idxTouch(kk) be a consecutive indices away from
        % idxPosStart 
        %   i.e., if idxPos = 44, 
        %         then idxTouch(kk) = 44,45,46,47,...
        %                   when kk =  1, 2, 3, 4,...
        % 
        %         idxTouch(kk) = 47 when kk = 2 means that idxTouch has skipped
        %         45 and 46. Therefore the while-loop breaks
        kStart = idxDiff(diffTemp == nanmax(diffTemp(diffTemp<=0)));%same as contactWidthIdxStop
        kk = kStart;
        while((kk <= numel(idxTouch)) && (idxTouch(kk) == (contactWidthIdxStart+(kk-kStart))))
            contactWidthIdx(2) = idxTouch(kk);%locates the single peak in profileGreen this farthest to the RIGHT of the contactPt at or around contactWidthIdxStart
            kk = kk + 1;
        end

        kk = kStart;
        while((kk >= 1) && (idxTouch(kk) == (contactWidthIdxStart-(kk-kStart))))
            contactWidthIdx(1) = idxTouch(kk);%locates the single peak in profileGreen this farthest to the LEFT of the contactPt at or around contactWidthIdxStart
            kk = kk - 1;
        end

    else %two indices closest to tempXY/XZ/YZ(1,:) == 0 but they're NOT the same index. Therefore, contactPt is NOT located around tempXY/XZ/YZ(1,:) == 0 
        %contactWidthIdxStart
        %           |
        %-----(-)<--|-----|------>(+)---|--
        %                 0             |
        %                       contactWidthIdxStop

        boolContactPt = 0;%NO, contactPt is NOT located around tempXY/XZ/YZ(1,:) == 0
        contactWidthIdx = [contactWidthIdxStart,contactWidthIdxStop];
    end
end


if((numel(diffTemp(diffTemp<=0)) > 0)  && (numel(diffTemp(diffTemp>=0)) == 0))%the ONLY single peaks in profileGreen, where profileGreen > profileRed, are located to the LEFT of contactPt
    idxDiff = 1:numel(diffTemp);
    contactWidthIdxStart = idxTouch(diffTemp == nanmax(diffTemp(diffTemp<=0)));
    if(nanmax(diffTemp(diffTemp<=0)) == 0) %diffTemp == distance relative to the contactPt located at tempXY/XZ/YZ(1,:) == 0
        %       contactWidthIdxStart
        %                 |
        % -----(-)<-------|------>(+)-----
        %                 0     

        boolContactPt = 1;%YES, contactPt is located around tempXY/XZ/YZ(1,:) == 0
    else
        %  contactWidthIdxStart
        %             |    
        % -----(-)<---|----|------>(+)-----
        %                  0     
        boolContactPt = 0;%NO, contactPt is NOT located around tempXY/XZ/YZ(1,:) == 0
    end

    %-----------(same 2 while-loops)----------
    kStart = idxDiff(diffTemp == nanmax(diffTemp(diffTemp<=0)));%same as contactWidthIdxStop
    kk = kStart;    
    while((kk <= numel(idxTouch)) && (idxTouch(kk) == (contactWidthIdxStart+(kk-kStart))))
        contactWidthIdx(2) = idxTouch(kk);%locates the single peak in profileGreen this farthest to the RIGHT of the contactPt at or around contactWidthIdxStart
        kk = kk + 1;
    end

    kk = kStart;
    while((kk >= 1) && (idxTouch(kk) == (contactWidthIdxStart-(kk-kStart))))
        contactWidthIdx(1) = idxTouch(kk);%locates the single peak in profileGreen this farthest to the LEFT of the contactPt at or around contactWidthIdxStart
        kk = kk - 1;
    end
end

if((numel(diffTemp(diffTemp<=0)) == 0)  && (numel(diffTemp(diffTemp>=0)) > 0))%the ONLY single peaks in profileGreen, where profileGreen > profileRed, are located to the RIGHT of contactPt
    idxDiff = 1:numel(diffTemp);
    contactWidthIdxStart = idxTouch(diffTemp == nanmin(diffTemp(diffTemp>=0)));
    if(nanmin(diffTemp(diffTemp>=0)) == 0) %diffTemp == distance relative to the contactPt located at tempXY/XZ/YZ(1,:) == 0
        %       contactWidthIdxStart
        %                 |
        % -----(-)<-------|------>(+)-----
        %                 0     
        boolContactPt = 1;%YES, contactPt is located around tempXY/XZ/YZ(1,:) == 0
    else
        %              contactWidthIdxStart
        %                       |
        % -----(-)<-------|-----|->(+)-----
        %                 0     

        boolContactPt = 0;%NO, contactPt is NOT located around tempXY/XZ/YZ(1,:) == 0
    end

    %-----------(same 2 while-loops)----------
    kStart = idxDiff(diffTemp == nanmin(diffTemp(diffTemp>=0)));%same as contactWidthIdxStop
    kk = kStart;    
    while((kk <= numel(idxTouch)) && (idxTouch(kk) == (contactWidthIdxStart+(kk-kStart))))
        contactWidthIdx(2) = idxTouch(kk);%locates the single peak in profileGreen this farthest to the RIGHT of the contactPt at or around contactWidthIdxStart
        kk = kk + 1;
    end

    kk = kStart;
    while((kk >= 1) && (idxTouch(kk) == (contactWidthIdxStart-(kk-kStart))))
        contactWidthIdx(1) = idxTouch(kk);%locates the single peak in profileGreen this farthest to the LEFT of the contactPt at or around contactWidthIdxStart
        kk = kk - 1;
    end
end
end