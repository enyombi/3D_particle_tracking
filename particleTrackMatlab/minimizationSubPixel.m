function [dp, RsquaredMin] =  minimizationSubPixel(centerIdx,centerDp,K,D,w,offset,cutoutSect,stepSize,deltaRcutoff,sampleSize)
% Author: Eru K.
% Data: 24-July-2014

% Objective: to return the sub-pixel displacement in row-, column-,
% slice-direction given by 'dp' from the sub-pixel position 'centerDp'

% 'dp' is located within a radial distance 'deltaRcutoff' from 'centerDp'
% BUT the position (i.e., row, col, slice) of 'dp' in the voxel is given
% using a coordinate sytem that is FIXED at 'center', i.e., [0,0,0]. So
% 'dp' is ALWAYS displacement wrt to the 'center' 
  
% RsqauredMin = sum of squares at dp

% Visual...
%                                   _______________________    
%                    .             / |                    /|
%                   .             /--|-------------------/ | 
%            (-1,-1,-1+stepSize) /:  |      :           /  |
%                               /______________________/   |
%                    (-1,-1,-1) | :  |      :          |   |   
%                               | :  |      :     c    |   |
%          (-1+stepSize,-1,-1)  | :  |      :          |   |
%                    .          | :  |      x----------|-- |
%                    .          | :  |______:__________|___|
%                    .          | : /       :          |   /(1,1,1)
%                               |  /--------:          |  /
%                               | /        /           | /
%                               |_________/____________|/
% 
% note: (-1,-1,-1) and (1,1,1) are just position in the cartesian
% coordinate system overlaid on the voxel
 
% center = ind2sub(size(cutoutSect),centerIdx);
% center = center of the particle AT PIXEL resolution in 'cutoutSect'
%      note: 'center' is located at the center of the voxel 
%      (i.e, pt 'x' where displacement = [0,0,0]
% 
% centerDp = subPixel components for the center of the particle 
%      'centerDp' is pt c or anywhere inside voxel centered around pt x
% 
% stepSize = spacing between each sub-pixel unit in the voxel  
%       note: stepSize basic defines the tics or spacing for a Cartesian 
%       coordinate system centered on pt x 

spacing = -1:stepSize:1; %this is x,y,z-tic marks for the Cartesian Coord system of the voxel, 'center' = [0,0,0] so dp is ALWAYS displacement from 'center'

[DeltaR,centerDpPos] = peakPlacementSubPixel(centerDp,stepSize); %DeltaR = matrix of radial distances from 'centerDp'; 
% note: the subscripted location (i.e., row, column, slice) 'centerDp' in
% DeltaR is given by 'centerDpPos'
% 
% DeltaR(centerDpPos(1),centerDpPos(2),centerDpPos(3)) = 0 
% 
% centerDp = [spacing(centerDpPos(1)),spacing(centerDpPos(2)),spacing(centerDpPos(3))]
%          = sub-pixel displacment from 'center' of voxel at (0,0,0)
% 
% Clarification: 
%   centerDp = displacement in row-,column-, and slice-direction from
%              'center' 
% 
%   centerDpPos = subscripted location of the 'centerDp' in 'DeltaR' 
%                 
% Important Note: DeltaR is NOT at the center of the voxel (i.e, pt x or
% 'center') but at some other position in SAME voxel pictured above.  This
% position is given by 'centerDpPos'. So, 
% 
%   DeltaR(centerDpPos(1),centerDpPos(2),centerDpPos(3)) =  0
% 
% because of roundoff error you can NOT compute the subscripted position of
% 'centerDp' using find()...
% 
%       centerDpPos(1) = find(spacing == centerDp(1)) 
%       centerDpPos(2) = find(spacing == centerDp(2)) 
%       centerDpPos(3) = find(spacing == centerDp(3)) 
% 
% ...instead you must peakPlacementSubPixel().m
bkgdVal = max(size(cutoutSect))*(10^9); %choose a background value well above any max Rsquared value so that Rsquared values are easily identified (note: after 1000 iterations is Rsquared ~10^5) 
Rsquared = bkgdVal*ones(size(DeltaR),'single');

% first compute Rsquared at the current displacement 'centerDp' whose
% subscripted position is given by 'centerDpPos'... 
[voronoiVol, Rmask] = peakPlacement(centerIdx,size(cutoutSect),D,w,centerDp);
mask = getCalcImg(Rmask,K,D,w,offset,voronoiVol);
sphericalRegion = find(voronoiVol == 1);
error = (cutoutSect(sphericalRegion)-mask(sphericalRegion)).^2;
Rsquared(centerDpPos(1),centerDpPos(2),centerDpPos(3)) = sum(error(:));
disp(sprintf('\t(%d of %d) Rsquared = %d, dp = [%.04f,%.04f,%.04f],  DeltaR = %.05f',0,sampleSize,Rsquared(centerDpPos(1),centerDpPos(2),centerDpPos(3)),centerDpPos(1),centerDpPos(2),centerDpPos(3),DeltaR(centerDpPos(1),centerDpPos(2),centerDpPos(3))));

% ...then compute Rsquared at a number (given by 'sampleSize') of randomly 
% selected postions less than a radial distance 'deltaRcutoff' from
% 'centerDp'. Because 'sampleRegion' has already been randomized just pick
% the first few elements in 'sampleRegion' 
%   i.e., pos == sampleRegion(1:sampleSize))
sampleRegion = find(DeltaR(:) < deltaRcutoff);
sampleRegion = sampleRegion(randperm(length(sampleRegion))); %randomize arrangement of indices in 'sphericalRegion' so that random sampling from 'sphericalRegion' is really RANDOM

[pos(:,1), pos(:,2), pos(:,3)] = ind2sub(size(DeltaR),sampleRegion(1:sampleSize)); %pos == subscripted position(i.e. pos(:,1) = row, pos(:,2) = col, pos(:,3) = slice)
for num = 1:sampleSize
    [voronoiVol, Rmask] = peakPlacement(centerIdx,size(cutoutSect),D,w,[spacing(pos(num,1)),spacing(pos(num,2)),spacing(pos(num,3))]);
    mask = getCalcImg(Rmask,K,D,w,offset,voronoiVol);
    sphericalRegion = find(voronoiVol == 1);
    error = (cutoutSect(sphericalRegion)-mask(sphericalRegion)).^2;
    Rsquared(pos(num,1),pos(num,2),pos(num,3)) = sum(error(:));
    disp(sprintf('\t(%d of %d) Rsquared = %d, dp = [%.04f,%.04f,%.04f],  DeltaR = %.05f',num,sampleSize,Rsquared(pos(num,1),pos(num,2),pos(num,3)),spacing(pos(num,1)),spacing(pos(num,2)),spacing(pos(num,3)),DeltaR(pos(num,1),pos(num,2),pos(num,3))));
end %20.4849seconds
RsquaredIdx = find(Rsquared(:) == min(Rsquared(:))); %the position of the min. val. should also be the position of the center of a particle
RsquaredMin = Rsquared(min(RsquaredIdx)); %min(RsquaredIdx) in case there are multiple minimum values of Rsquared.

[rowDp,colDp,sliceDp] = ind2sub(size(Rsquared),RsquaredIdx);
dp = [spacing(rowDp),spacing(colDp),spacing(sliceDp)];

disp(sprintf('RsquaredMin = %d, dp = [%.04f,%.04f,%.04f]',RsquaredMin,dp(1),dp(2),dp(3)));
end