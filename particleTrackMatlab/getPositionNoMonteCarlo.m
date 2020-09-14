function [peakNew, RsquaredNew] = getPositionNoMonteCarlo(topLeftCorner,center,K,D,w,offset,cutoutSect)
% Author: Eru K.
% Date: 14-July-2014
% adopted from getPosition.m
% Objective: to get the position of a particle by minimizing the error
% (really the sum of errors squared, Rsquared) between the mask image of   
% a particle given by ssf().m and the real image of particle in 'Inorm' 

% %    2________3
% %    /.       /|
% %  1/_.______/ |
% %   | .      4 |
% %   | .  x   | |
% %   | 6......|./7
% %   |._______|/
% %   5         8  
% % 
% % x = position of particle input, which is located at the center of 
% %     the cutoutSect

%minimizationPixelRes.m returns the position, 'center', of the lowest
%Rsquared value computed from ALL pixel positions in 'regionSize' 

%regionSize == radius of the spherical region surrounding 'center' in \
% 'cutoutSect' from which you're sampling
% note:  D = 85.6162  computes Rsquared for each of pixel position of a
% spherical region that is roughly 9x9x9pixels or radius of 3pixels 
% around 'center' 

regionSize = D*0.05; 
[centerNew, RsquaredNew] = minimizationPixelRes(center,K,D,w,offset,cutoutSect,regionSize); 
peakNew = topLeftCorner + (centerNew - ones(size(topLeftCorner),'single')); %the position of 'center' scaled to fit the fully-resolved image, 'Inorm'

% center(min(idxMin),:) is an actual location. But you need
% the displacement or distance between center(min(idxMin),:) and the 
% top left Corner of the cutoutSect to locate center(min(idxMin),:) in the 
% larger image from which cutoutSect is taken. To rescale
% center(min(idxMin),:) to the larger image.

%Visual: 
% center is located at pt 'x' 
% topLeftCorner_of_cutoutSect = [1,1,1] ALWAYS!
% x = topLeftCorner + distanceToCenter;
%   where, distanceToCenter = center - [1,1,1] because MATLAB starts indexing/counting at 1 
% 
%          _________
%         /.       /|
% [1,1,1]/_.______/ |
%        | .      | |
%        | .  x   | |
%        | .......|./
%        |._______|/
% 
% e.g., 
% Suppose topLeftCorner = [1,1,1] and center = [3,3,3]
%   distanceToCenter = center - [1,1,1]
%   distacneToCenter = [3,3,3] - [1,1,1] = [2,2,2]
% 
%   x = topLeftCorner_of_cutoutSect + distanceToCenter
%     = [1,1,1] + [2,2,2]
%     = [3,3,3]
% 
% Q.E.D., x = center

% So now in the larger image... 
% 
%        /
%       /      
%[1,1,1]_____________________________
%      |                      
%   topLeftCorner/_____/  |
%      |         |     |  |
%      |         |  x  | /
%      |         |_____|/
%      |
%      |
% 
% peakNew = topLeftCorner + distToCenter
% if topLeftCorner = [124 423 111] then,
% peakNew = [124,423,111] + [3,3,3]
%         = [127,426,114]
end