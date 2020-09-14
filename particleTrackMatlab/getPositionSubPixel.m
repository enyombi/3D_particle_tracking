function [peakNew, Rsquared] = getPositionSubPixel(topLeftCorner,centerIdx,K,D,w,offset,cutoutSect,rangeColormap)
% Author: Eru K.
% Data: 24-July-2014
% adopted from testSubPixelRes3.m
% Objective: find the center of a particle at sub-pixel resolution

stepSize = single(1/rangeColormap);%max. resolution

iterationsM = 11; %iterations of the Monte Carlo optimization rountine in minimizationSubPixel().m
Rsquared = zeros([iterationsM, 1],'single');
dp = zeros([iterationsM, 3],'single'); %sub-pixel displacement, note: dp(1,:) = [0,0,0]
deltaRcutoff = linspace(1,10*stepSize,iterationsM); %in fractions of pixels

% compute Rsquared at dp = [0,0,0]
[voronoiVol, Rmask] = peakPlacement(centerIdx,size(cutoutSect),D,w,dp(1,:)); 
mask = getCalcImg(Rmask,K,D,w,offset,voronoiVol);
sphericalRegion = find(voronoiVol == 1);
error = (cutoutSect(sphericalRegion)-mask(sphericalRegion)).^2;
Rsquared(1) = sum(error(:));

sampleSize = 100; %number of time to sample from the population given by sampleRegion

for numM = 2:iterationsM
    disp(sprintf('\n---------------(START, iteration %d)---------------',numM));
    [dp(numM,:), Rsquared(numM)] = minimizationSubPixel(centerIdx,dp(numM-1,:),K,D,w,offset,cutoutSect,stepSize,deltaRcutoff(numM-1),sampleSize);
    disp(sprintf('---------------(END, iteration %d)---------------',numM));
end

idxMin = find(Rsquared(:) == min(Rsquared(:)));

[center(1,1),center(1,2),center(1,3)] = ind2sub(size(cutoutSect),centerIdx);
centerNew = center + dp(min(idxMin),:); %particle position in 'cutoutSect at sub-pixel resolution
peakNew = (topLeftCorner + centerNew) - ones(size(topLeftCorner),'single'); %the position of 'centerNew' scaled to fit the fully-resolved image, 'Inorm'
end