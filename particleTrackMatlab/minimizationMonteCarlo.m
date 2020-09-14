function [posRsquaredMin, RsquaredMin] = minimizationMonteCarlo(center,K,D,w,offset,cutoutSect,regionSize,sampleSize)
% Author: Eru K.
% Date: 30-Jun-2014
% objective: To compute Rsquared, i.e. the sum of residuals squared, for a
% random sampling of pixel-positions that are a distance D/2 from 'center'.
% And return the lowest value of Rsquared and its position.

% note: the residuals are the differences between the mask image (ssf().m)
% located at any given pixel that is a distance D/2 from 'center' and the
% real image 'cutoutSect'

% center = position of the particle center subscripted i.e., center = row,
% col, slice]
% D = particle diameter
% w = thickness of particle boundary
% cutoutSect = section cutout from 'Inorm'
% regionSize = size of the spherical region surrounding 'center' defined by
%              a distance from 'center' 
% sampleSize = number of samples to be collected from 'regionSize'
% (i.e.,Monte Carlo sample size)

bkgdVal = max(size(cutoutSect))*(10^9); %choose a background value well above the max Rsquared value so that Rsquared values are easily identified (note: after 1000 iterations is Rsquared ~10^5) 
Rsquared = bkgdVal*ones(size(cutoutSect),'single');

% first compute Rsquared at the current position....
centerIdx = sub2ind(size(cutoutSect),center(1),center(2),center(3));
[voronoiVol,R] = peakPlacement(centerIdx,size(cutoutSect),D,w); %'peakPlacement' returns the spherical region centered on the index located at sphericalRegion(np) for drawing the particle. 'voronoiVol' is the spherical region and 'R' is the radial values to be used in ssf().m
calcCutoutImg = getCalcImg(R,K,D,w,offset,voronoiVol);
newSphericalRegion = find(voronoiVol == 1);
error = (cutoutSect(newSphericalRegion)-calcCutoutImg(newSphericalRegion)).^2;
Rsquared(centerIdx) = sum(error(:)); 

% now compute Rsquared at random positions at distance less than or equal
% to D/2 from the particle's center
[row, col, slice] = ndgrid(1:size(cutoutSect,1),1:size(cutoutSect,2),1:size(cutoutSect,3));

rowTemp = row - center(1); %y-displacement from center(np, )
colTemp = col - center(2); %x-displacement from center(np, )
sliceTemp = slice - center(3); %z-displacment from center(np, )

rTemp = sqrt(rowTemp.^2 + colTemp.^2 + sliceTemp.^2); %Eqn of Sphere: (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = Radius^2, where center = (x0,y0,z0)

sphericalRegion = find(rTemp(:) <= regionSize); % 'sphericalRegion' == sampleRegion; identifies the indices of elements in 'cutoutSect' contained in a spherical region surrounding 'peak'.  
sphericalRegion = sphericalRegion(randperm(length(sphericalRegion))); %randomize arrangement of indices in 'sphericalRegion' so that random sampling from 'sphericalRegion' is really RANDOM

idx = randi(length(sphericalRegion),[sampleSize 1]); % Randomly select round((length(sphericalRegion)*fract)) locations in 'sphericalRegion'

disp(sprintf('Rsquared = %d for index %d out of %d',Rsquared(centerIdx),0,sampleSize))
for np = 1:length(idx)
    [voronoiVol, R] = peakPlacement(sphericalRegion(np),size(cutoutSect),D,w); %'peakPlacement' returns the spherical region centered on the index located at sphericalRegion(np) for drawing the particle. 'voronoiVol' is the spherical region and 'R' is the radial values to be used in ssf().m
    calcCutoutImg = getCalcImg(R,K,D,w,offset,voronoiVol);
    newSphericalRegion = find(voronoiVol == 1);
    error = (cutoutSect(newSphericalRegion)-calcCutoutImg(newSphericalRegion)).^2;
    Rsquared(sphericalRegion(np)) = sum(error(:));
    disp(sprintf('Rsquared = %d for index = %d; iter num = %d',Rsquared(sphericalRegion(np)),sphericalRegion(np),np))
end

RsquaredIdx = find(Rsquared(:) == min(Rsquared(:))); %the position of the min. val. should also be the position of the center of a particle
[posRsquaredMin(1), posRsquaredMin(2), posRsquaredMin(3)] = ind2sub(size(cutoutSect),RsquaredIdx);% posRsquaredMin = center;
RsquaredMin = Rsquared(min(RsquaredIdx)); %min(RsquaredIdx) in case there are multiple minimum values of Rsquared.
disp(sprintf('RsquaredMin = %d',RsquaredMin))
end

% % -----------testing....
% K = Kactual;
% D = Dactual;
% w = wActual;
% a = find(validPks(:,3) == 619);
% peak = validPks(a,:);
% 
% dist = 2*D;
% % vertices of cutoutSect
% v(1,:) = peak - dist; %topFrontCorner
% v(7,:) = peak + dist; %bottomBackCorner
% 
% v(2,:) = [peak(1)-dist, peak(2)-dist, peak(3)+dist];
% v(3,:) = [peak(1)-dist, peak(2)+dist, peak(3)+dist];
% v(4,:) = [peak(1)-dist, peak(2)+dist, peak(3)-dist];
% 
% v(5,:) = [peak(1)+dist, peak(2)-dist, peak(3)-dist];
% v(6,:) = [peak(1)+dist, peak(2)-dist, peak(3)+dist];
% v(8,:) = [peak(1)+dist, peak(2)+dist, peak(3)-dist];
% 
% % now make certain that the rectangular box around the peak does not
% % exceed the bounds of the 3D image...
% for np = 1:length(v) %particle #
%     for dim = 1:3 %each particle has 3 dimensions/coords
%         if(v(np,dim) <= 0)
%             v(np,dim) = 1;
%         end
% 
%         if(v(np,dim) > size(Inorm,dim))
%             v(np,dim) = size(Inorm,dim);
%         end
%     end
% end
% 
% v = round(v);
% cutoutSect = Inorm(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
% center = ((size(cutoutSect)-1)/2)+1;
% regionSize = D;
% sampleSize = 10;