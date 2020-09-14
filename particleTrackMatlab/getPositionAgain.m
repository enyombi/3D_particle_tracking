function [peakNew, Rsquared] = getPositionAgain(topLeftCorner,centerCutoutSect,K,D,w,offset,cutoutSect)
% Author: Eru K.
% Date: 29-July-2014
% Objective: to get the position of a particle by minimizing the error
% (really the sum of errors squared, Rsquared) between the mask image of   
% a particle given by ssf().m and the real image of particle in 'Inorm' 

% This is exactly the same at getPosition.m except it conducts the Monte
% Carlo Minimization over a smaller region using Inorm instead of Inorm2,
% Inorm3, etc. 

% The idea is that peaks located by getPosition.m are on a particle
% but because getPosition.m located these peaks on an a version of
% 'Inorm'...
% 
% (e.g., Inorm2, Inorm3, etc.  
%       Inorm2 has particle located in the first iterations of particle
%       tracking on Inorm removed 
%   
%       Inorm3 has particles located in the first and 2nd iterations of the
%       particle tracking algorithm removed)
% 
% ...with some particles removed, the position of the peak may not be
% located precisely at the particle's center. Therefore, run
% getPositionAgain.m using Inorm on a smaller cutoutSect and using smaller
% regionSize so that the peak is moved on the SAME particle to locate the
% center better.

% You know that the peak (i.e., 'centerCutoutSect') is on the particle and
% very close to that particle center you just want to move it and get a
% more accurate estimate of the particle center using 'Inorm' instead of
% some subsequent iteration of Inorm that has particle removed (e.g.,
% Inorm2, Inorm3, etc.)

% -----------------Monte Carlo optimization-----------------
%     (i.e., minimize Rsquared to find particle position)

% Each iteration computes Rsquared at a number of pixel positions that are 
% a distance of 'regionSize' from 'center'. ('sampleSize' = the number of 
% pixel positions randomly chosen from all pixels within a distance 
% 'regionSize' from 'center')
% 
% The  smallest Rsquared value and its position are returned by
% getPixelRes().m
% 
% Note: the position of the smallest Rsquared value is the particle's 
% 'center'
iterationsM = 3; %iterations of the Monte Carlo optimization rountine in getPixelRes().m
regionSize = linspace(D*0.07,D*0.05,iterationsM);%make regionSize fairly small to make sure that you stay on the same particle 

Rsquared = zeros([iterationsM, 1],'single');
sampleSize = 100; %it works for sampleSize = 100 and iterationsM = 10

% regionSize == radius of the spherical region surrounding the pixel
% position  (i.e., 'center') in 'cutoutSect' from which you're sampling. 
% Each iteration of of the Monte Carlo optimization (i.e., getPixelRes().m)
% returns the position of the smallest Rsquared value and brings you closer
% to the particle's center. To enhance efficiency, the size of the sampling
% region (i.e., regionSize) also decreases each iteration until reaching
% the ideal size of a particle, D/2. The idea is that reducing 'regionSize'
% with each iteration causes the entire Monte Carlo optimization routine to 
% converge on correct particle position. 
[center(1,:), Rsquared(1)] = minimizationMonteCarlo(centerCutoutSect,K,D,w,offset,cutoutSect,regionSize(1),sampleSize); 
for num = 2:iterationsM
    [center(num,:), Rsquared(num)] = minimizationMonteCarlo(center(num-1,:),K,D,w,offset,cutoutSect,regionSize(num-1),sampleSize);
end

% Monte Carlo Minimization should put you very close to the 'true' particle
% position at pixel resolution. Now compute Rsquared for ALL locations 
% VERY close to and w/in a distance of 'regionSize' from center(idxMin,:).
% Retain the position with the lowest Rsquared value so
% center(iterationsM+1,:) should be the location with the smallest Rsquared
% value
idxMin = find(Rsquared(:) == min(Rsquared(:)));
regionSize = D*0.03; %tests a box of 9x9x9pixels around center(min(idxMin),:) for D = 85.6162
[centerNew, RsquaredNew] = minimizationPixelRes(center(min(idxMin),:),K,D,w,offset,cutoutSect,regionSize); 

center = vertcat(center,centerNew);
Rsquared = vertcat(Rsquared,RsquaredNew);

idxMin = find(Rsquared(:) == min(Rsquared(:)));
peakNew = topLeftCorner + (center(min(idxMin),:) - ones(size(topLeftCorner),'single')); %the position of 'center' scaled to fit the fully-resolved image, 'Inorm'

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