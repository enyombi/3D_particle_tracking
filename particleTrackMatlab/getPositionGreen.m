function [centerNew, radiusNew, validPtsEgdes] = getPositionGreen(topLeftCorner,center,cutoutSectGreen,cutoutSectRed,particleRadius)
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

% pts = single(pts);

% the radial distance between points on opposite sides/poles of the
% sphere should be the same even if the sphere is ellipsoidal and not
% spherical


% make sure the vertices of the cube circumscirbing the sphere are inside 'cutoutSect*'
for npp = 1:size(pts,1)
    for jj = 1:3
        if(pts(npp,jj) < 1)
            pts(npp,jj) = 1;
        end
        
        if(pts(npp,jj) > size(cutoutSectGreen,jj))
            pts(npp,jj) = size(cutoutSectGreen,jj); 
        end
    end     
end


radiusApprox = nan(size(pts,1),1);
for npp = 1:size(pts,1)
    [profileRedTemp, diagPlaneRed] = getLineProfile(center,pts(npp,:),cutoutSectRed);
    [profileGreenTemp, diagPlaneGreen] = getLineProfile(center,pts(npp,:),cutoutSectGreen);
    
    localMaximaRedIdx = getLocalMaxima(profileRedTemp); %location of the peaks in profileGreenTemp
    localMaximaGreenIdx = getLocalMaxima(profileGreenTemp); %location of the peaks in profileGreenTemp
    if(~isempty(localMaximaGreenIdx))
        %assume closest peak in 'profileGreenTemp' to 'center' is the edges
        %of the particle
        radiusApprox(npp) = min(localMaximaGreenIdx);
        temp(npp) = min(localMaximaGreenIdx);
        
        %Now correct your assumption!  By assuming that the largest peak in
        %'profileRedTemp' is just outside the fluorescent green edges of
        %the particle.  Therefore, the first green peak just before the max
        %peak in 'profileRedTemp' is the actual edge of the particle 
        
        %Peaks in 'profileRedTemp' are NOT trustworthy because they
        %measure the red fluorescence surrounding the particle (i.e.,
        %fluorescence that is exlcuded from the space that the particle
        %occupies). Noise, artefacts, point spread funct, etc. may show red
        %fluorscence inside the particle which is not accurate!
        %Therefore, only proceed with this 2nd assumption if the largest
        %peak in 'profileRedTemp' is located at a radial distance greater
        %than the smallest peak in 'profileGreenTemp'  
        if(max(localMaximaRedIdx) > min(localMaximaGreenIdx)) 
            tempProfileRedMaxima = profileRedTemp(localMaximaRedIdx); %max pixel values in profileRedTemp
            largestRed = find(tempProfileRedMaxima == max(tempProfileRedMaxima));
            largestMaxRedIdx = localMaximaRedIdx(largestRed);

            %modify localMaximaGreenIdx to only list those peaks in
            %profileGreenTemp that appear before the max peak in profileRedTemp
            localMaximaGreenIdxNew = localMaximaGreenIdx;
            localMaximaGreenIdxNew(localMaximaGreenIdx > largestMaxRedIdx) = [];
            
            if(~isempty(localMaximaGreenIdxNew))
                tempProfileGreenMaxima = profileGreenTemp(localMaximaGreenIdxNew); %max pixel values in profileGreenTemp that appear before the max pixel value in profileRedTemp
                largestGreen = find(tempProfileGreenMaxima == max(tempProfileGreenMaxima));
                largestMaxGreenIdx = localMaximaGreenIdxNew(largestGreen);
                
                %If the largest peak in 'profileGreenTemp' just before the
                %largest peak in 'profileRedTemp' is nearly the same
                %intensity as the peak in 'profileGreenTemp' that is
                %closest to 'center' then first peak is likely to be the
                %edge of the particle around 'center' and the other peak is
                %just the edge of a neighboring particle. 
                %So only replace radiusApprox if
                %'profileGreenTemp(largestMaxGreenIdx)' is >
                %0.3*profileGreenTemp(min(largestMaxGreenIdx))
                if((profileGreenTemp(largestMaxGreenIdx)-profileGreenTemp(min(localMaximaGreenIdx)))/profileGreenTemp(min(localMaximaGreenIdx)) > 0.3)
                    radiusApprox(npp) = largestMaxGreenIdx;
                end
            end
        end
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


% replace ptsEdges with NaN if they are outside the image
for npp = 1:size(ptsEdges,1)
    for jj = 1:3
        if(ptsEdges(npp,jj) < 1)
            ptsEdges(npp,:) = nan;
            radiusApprox(npp) = nan; %the edge of the particle is outside of the image so disregard radiusApprox(npp) as well
        end
        
        if(ptsEdges(npp,jj) > size(cutoutSectGreen,jj))
            ptsEdges(npp,:) = nan; 
            radiusApprox(npp) = nan; %the edge of the particle is outside of the image so disregard radiusApprox(npp) as well
        end
    end     
end

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
correctionFactor = round(0.75*nanmean(radiusApprox(:)));

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
possibleIdx = sub2ind(size(cutoutSectGreen),rowRange(1),colRange(1),sliceRange(1));
for row = rowRange(1):rowRange(2)
    for col = colRange(1):colRange(2)
        for slice = sliceRange(1):sliceRange(2)
            possibleIdx = vertcat(possibleIdx,sub2ind(size(cutoutSectGreen),row,col,slice));%vertcat() takes ~25sec BUT union() takes ~5min
        end
    end
end
possibleIdx(1) = [];

idxSample = randperm(numel(possibleIdx),numTry);
% 1st guess is 'center' and 'radiusApprox':
centerTest(1,:) = center;
radius(:,1) = radiusApprox;
rangeRadii(1) = nanmax(radiusApprox(:))-nanmin(radiusApprox(:));
for num = 2:numel(idxSample)        
    [rowTemp,colTemp,sliceTemp] = ind2sub(size(cutoutSectGreen),possibleIdx(idxSample(num)));
    centerTest(num,:) = [rowTemp,colTemp,sliceTemp];
    
    for npp = 1:size(ptsEdges,1)
        radius(npp,num) = sqrt(((rowTemp - ptsEdges(npp,1)).^2) + ((colTemp - ptsEdges(npp,2)).^2) + ((sliceTemp - ptsEdges(npp,3)).^2));
    end
    
    rangeRadii(num) = nanmax(radius(:,num))-nanmin(radius(:,num));
end

for num = 1:numTry
    avg(num) = nanmean(radius(:,num));
    standDev(num) = nanstd(radius(:,num),1);
end

npBest = find(rangeRadii == nanmin(rangeRadii(:)));

centerNew = topLeftCorner + centerTest(npBest,:) - ones(size(center));
radiusNew = nanmean(radius(:,npBest)); %average of up to 12 radii measurements

temp = isnan(radius(:,npBest));
validPtsEgdes = numel(find(temp == 0));
end