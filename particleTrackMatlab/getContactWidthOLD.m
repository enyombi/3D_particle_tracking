function contactWidthIdx = getContactWidth(posR,contactRegion,boolMaxima)
% Author: Eru K.
% objective: To acquire the width of the contact region. 
%   (This is the number of locations adjacent to the point of contact where
%   there is a single peak in profileGreen)
% 
% posR == position wrt contact pt. 'x'
% 
% contactRegion(1,:) == position of peaks in InormGreen near particle np
% 
% contactRegion(2,:) == position of peaks in InormGreen near particle
%                          npp
% 
% boolMaxima(jj,1) == 1 if peak in InormGreen measured near particle np is
%                     real
%                  == 0 if peak is NOT real (so, that the edge of particle
%                     np is fake
% 
% boolMaxima(jj,2) <--same as boolMaxima(jj,1) just for particle npp

% ref: Landau and Lifshitz, pgs. 31-34
% 
% particles np and npp
% 'x' = presumed center of the contact region
% w = contactWidth
% 
% sideview:
%             particle npp,center
%                     A
%                     |
%         \           |           /<---contactRegion(2,:)
%          \          |          /
%           \_________x_________/
%           /|        |        |\
%          / |        |        | \
%         /  |        |        |  \<---contactRegion(1,:) 
%            |        |        |
%            |        |        |
%            | (-)<---0--->(+) |
%                     |
%               find(posR==0)
%                     |
%                     |
%                     V
%               particle np,center
% 
% note: frontview of the contact Region is a circle or an ellipse
contactRegion(boolMaxima(:,1)==0,1) = nan; %replace real-values with nan-values for intensity in profileGreen that are NOT peaks and therefore do NOT represent the edges of a particle
contactRegion(boolMaxima(:,2)==0,2) = nan;
coincidence = contactRegion(:,1)-contactRegion(:,2);%co-incident peaks in profileGreen for np and npp identify the width of the contactRegion

idx = 1:numel(coincidence); %all indices
idx(boolMaxima(:,1)==0) = nan;
idx(boolMaxima(:,2)==0) = nan;
idx(isnan(idx)==1) = [];

%IMPORTANT: a co-incident btwn contactRegion(:,1) and
%contactRegion(:,2) includes only those values of coincidence == 0
%   valCoincidence = coincidence(idx);
%   idxCoincidence = idx(valCoincidence == 0);

idxCoincidence = idx(coincidence(idx) == 0); %all co-indicent indices
diffTemp = abs(posR(coincidence(idx) == 0)); %distance each index is from posR == 0. Equivalent to abs(idxCoincidence-find(posR==0));

contactWidthStart = idxCoincidence(diffTemp==nanmin(diffTemp(:)));

% note: diffTemp and idxCoincidence are the same length so the location of
% the min. value in diffTemp also identifies its corresponding index

% [idxCoincidence',diffTemp']
% 
% ans =
% 
%     29    10
%     30     9
%     33     6
%     36     3
%     38     1
%     39     0 <---contactWidthStart
%     41     2

% find the co-incidence farthest to the left/(-)-direction of the
% find(posR==0)
contactWidthStop = nan([1,2]);
kk = find(posR==0);
while(((abs(contactRegion(kk,1))<2) || (isnan(contactRegion(kk,1)) == 1)) && ((abs(contactRegion(kk,2))<2) || (isnan(contactRegion(kk,2)) == 1)) && (kk>1)) %proceed as long as abs(contactRegion(kk,1))<2 or nan (samme goes for contactRegion(kk,2)) AND 1<=kk<=find(posR==0) 
    if(coincidence(kk) ==0)
        contactWidthStop(1) = kk;
    end
    kk = kk - 1;
end

kk = find(posR==0);
while(((abs(contactRegion(kk,1))<2) || (isnan(contactRegion(kk,1)) == 1)) && ((abs(contactRegion(kk,2))<2) || (isnan(contactRegion(kk,2)) == 1)) && (kk<numel(coincidence))) %proceed as long as abs(contactRegion(kk,1))<2 or nan (samme goes for contactRegion(kk,2)) AND find(posR==0)<=kk<=numel(coincidence) 
    if(coincidence(kk) ==0)
        contactWidthStop(2) = kk;
    end
    kk = kk + 1;
end

% contactWidthStart is ALWAYS populated by 1 or more real-values. But these
% are NOT always leftmost(-) and rightmost(+) coincidences 
%   i.e., sometimes: 
%       contactWidthStop(1) == nan 
%           or 
%       contactWidthStop(2) == nan 
% When this happens replaces either index contactWidthStop(1) or
% contactWidthStop(2) with contactWidthStart...
if(isnan(contactWidthStop(1))==1) %no coincidences found to the left/(-)-direction of find(posXY==0)
    contactWidthStop(1) = nanmin(contactWidthStart(:));%nanmin(contactWidthStart(:)) == leftmost coincidence that is closest to find(posR == 0) and farthest from contactWidthStop(2)
end

if(isnan(contactWidthStop(2))==1) %no coincidences found to the right/(+)-direction of find(posXY==0)
    contactWidthStop(2) = nanmax(contactWidthStart(:));%nanmax(contactWidthStart(:)) == rightmost coincidence that is closest to find(posR == 0) and farthest from contactWidthStop(1)
end

contactWidthIdx = contactWidthStop;
end