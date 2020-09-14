function [centerNew, radiusNew, validPtsEgdes] = getPositionGreen(topLeftCorner,center,cutoutSect,particleRadius)
% Author: Eru K.
% Date: 5-may-2015

% Objective: locates the center using the green fluorescent image
%   1) first find the fluorescent edges of a particle, 'ptsEdges'
%   2) measure the distance btwn 'ptsEdges' and an arbitrary 'centerTest'
%   
%   the position 'centerTest' at which the radial distance between each
%   'ptsEdges' are equal is the most precise approximation of the
%   particle's center.  This position is returned as 'centerNew'. The
%   average distance between 'centerNew' and each 'ptsEdges' is returned as
%   'radiusNew'

% -------------------------------x-----------------------------------
% Step1: find the fluorescent edges of the particle, i.e., 'ptsEdges'

% Draw a box around center, 'c'
% 
%     slices(+)
%    /
%   /
%   ----->cols(+)
%   |
%   |
%   rows(+)
% 
%    5______________4
%    /:     /'     /|
%   /_:__ 9/_'____/ | 
% 1/:_:__ /:_8___/| |
% |11_: _| c_'_ _||12
% |/:/_ _|/:/____|/ |
% | :   7| /10   ||/ 2
% |/____ |/______|/ 
% 3              6

% vertices of diagonal lines through ptC
pts(1,:) = center - particleRadius;
pts(2,:) = center + particleRadius;

pts(3,:) = [center(1)+particleRadius, center(2)-particleRadius, center(3)-particleRadius];
pts(4,:) = [center(1)-particleRadius, center(2)+particleRadius, center(3)+particleRadius];

pts(5,:) = [center(1)-particleRadius, center(2)-particleRadius, center(3)+particleRadius];
pts(6,:) = [center(1)+particleRadius, center(2)+particleRadius, center(3)-particleRadius];

% vertices of vertical and horizontal lines through ptC
pts(7,:) = [center(1), center(2), center(3)-particleRadius];
pts(8,:) = [center(1), center(2), center(3)+particleRadius];

pts(9,:) = [center(1)-particleRadius, center(2), center(3)];
pts(10,:) = [center(1)+particleRadius, center(2), center(3)];

pts(11,:) = [center(1), center(2)-particleRadius, center(3)];
pts(12,:) = [center(1), center(2)+particleRadius, center(3)];

pts = single(pts);

% the radial distance between points on opposite sides/poles of the
% sphere should be the same even if the sphere is ellipsoidal and not
% spherical

radiusApprox = nan(size(pts,1),1);
for npp = 1:size(pts,1)
    [profileGreenTemp, diagPlaneGreen] = getLineProfile(center,pts(npp,:),cutoutSect);
    localMaximaIdx = getLocalMaxima(profileGreenTemp); %location of the peaks in profileGreenTemp
    if(~isempty(localMaximaIdx))
        radiusApprox(npp) = min(localMaximaIdx);%radial distance from 'center' to the nearest peak in 'profileGreenTemp'.  This is presumably the edge of the particle
    end
end

% ptsEdges == points that are actually on the the fluorescent edges of the
% particle
ptsEdges = nan(size(pts));

ptsEdges(1,:) = center - radiusApprox(1);
ptsEdges(2,:) = center + radiusApprox(2);

ptsEdges(3,:) = [center(1)+radiusApprox(3), center(2)-radiusApprox(3), center(3)-radiusApprox(3)];
ptsEdges(4,:) = [center(1)-radiusApprox(4), center(2)+radiusApprox(4), center(3)+radiusApprox(4)];

ptsEdges(5,:) = [center(1)-radiusApprox(5), center(2)-radiusApprox(5), center(3)+radiusApprox(5)];
ptsEdges(6,:) = [center(1)+radiusApprox(6), center(2)+radiusApprox(6), center(3)-radiusApprox(6)];

% vertices of vertical and horizontal lines through ptC
ptsEdges(7,:) = [center(1), center(2), center(3)-radiusApprox(7)];
ptsEdges(8,:) = [center(1), center(2), center(3)+radiusApprox(8)];

ptsEdges(9,:) = [center(1)-radiusApprox(9), center(2), center(3)];
ptsEdges(10,:) = [center(1)+radiusApprox(10), center(2), center(3)];

ptsEdges(11,:) = [center(1), center(2)-radiusApprox(11), center(3)];
ptsEdges(12,:) = [center(1), center(2)+radiusApprox(12), center(3)];


% replace all ptsEdges with NaN if radiusApprox(npp) = nan
for npp = 1:size(ptsEdges,1)
    if(isnan(radiusApprox(npp)))
        ptsEdges(npp,:) = nan(1,3);
    end
end


%     slices(+)
%    /
%   /
%   ----->cols(+)
%   |
%   |
%   rows(+)
%                       
%              9    8
%        5_____:___/____4
%        /     :  /    /|
%       /      : /    / | 
%     1/______ :/____/  |
%  11_| _ _ _ _c__ _ | _|_ _12
%     |       /:     |  |
%     |      / :     | / 2
%     |_____/__:_____|/ 
%     3    /   :     6
%         /   10 
%        7      

% try re-assigning the center of the sphere to any pt. inside the cube that
% circumscribes the sphere whose edges are 'ptsEdges'

rowRange = round([nanmin(ptsEdges(:,1)),nanmax(ptsEdges(:,1))]);
colRange = round([nanmin(ptsEdges(:,2)),nanmax(ptsEdges(:,2))]);
sliceRange = round([nanmin(ptsEdges(:,3)),nanmax(ptsEdges(:,3))]);

% make the cube smaller...
% 
% _____rowRange(2)=143
%   |
%   |           _____rowRange(2)-correctionFactor = 143-10 = 133
%   |             |  
%   |             |
%   |             |  
%   |           __|__rowRange(1)+correctionFactor = 63+10 = 73   
%   |
% __|__
%      rowRange(1)=63 
% 
correctionFactor = round(0.5*nanmean(radiusApprox(:)));

if((rowRange(2) - correctionFactor) > (rowRange(1) + correctionFactor + 10)) %ensures that rowRange(2) - rowRange(1) > 10 after adding the correctionFactor Before applying the correctionFactor
    rowRange(2) = rowRange(2) - correctionFactor;
    rowRange(1) = rowRange(1) + correctionFactor;
end

if((colRange(2) - correctionFactor) > (colRange(1) + correctionFactor + 10)) %ensures that colRange(2) - colRange(1) > 10 after adding the correctionFactor Before applying the correctionFactor
    colRange(2) = colRange(2) - correctionFactor;
    colRange(1) = colRange(1) + correctionFactor;
end

if((sliceRange(2) - correctionFactor) > (sliceRange(1) + correctionFactor + 10)) %ensures that sliceRange(2) - sliceRange(1) > 10 after adding the correctionFactor Before applying the correctionFactor
    sliceRange(2) = sliceRange(2) - correctionFactor;
    sliceRange(1) = sliceRange(1) + correctionFactor;
end

numTry = min(6000,round(0.05*(rowRange(2)-rowRange(1))*(colRange(2)-colRange(1))*(sliceRange(2)-sliceRange(1)))); %No more than 6000 guesses for the particle's center
radius = nan(size(ptsEdges,1),numTry);
centerTest = nan(numTry,3);

% stats for radii...
rangeRadii = nan(numTry,1);

% possibleIdx = all the pixel locations of the cube that  
% circumscribes the sphere whose edges are 'ptsEdges'  
possibleIdx = sub2ind(size(cutoutSect),rowRange(1),colRange(1),sliceRange(1));
for row = rowRange(1):rowRange(2)
    for col = colRange(1):colRange(2)
        for slice = sliceRange(1):sliceRange(2)
            possibleIdx = union(possibleIdx,sub2ind(size(cutoutSect),row,col,slice));
        end
    end
end


idxSample = randperm(numel(possibleIdx),numTry);

for num = 1:numel(idxSample)        
    [rowTemp,colTemp,sliceTemp] = ind2sub(size(cutoutSect),possibleIdx(idxSample(num)));
    centerTest(num,:) = [rowTemp,colTemp,sliceTemp];
    
    for npp = 1:size(ptsEdges,1)
        radius(npp,num) = sqrt(((rowTemp - ptsEdges(npp,1)).^2) + ((colTemp - ptsEdges(npp,2)).^2) + ((sliceTemp - ptsEdges(npp,3)).^2));
    end
    
    rangeRadii(num) = nanmax(radius(:,num))-nanmin(radius(:,num));
end

npBest = find(rangeRadii == nanmin(rangeRadii(:)));

centerNew = topLeftCorner + centerTest(npBest,:) - ones(size(center));
radiusNew = nanmean(radius(:,npBest)); %average of up to 12 radii measurements

temp = isnan(radius(:,npBest));
validPtsEgdes = numel(find(temp == 0));
end