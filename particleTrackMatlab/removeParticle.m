function [diffImg, Rsquared, peakNum] = removeParticle(peak,K,D,w,offset,cutoutSect,np)
% Author: Eru K.
% date: 8-July-2014

% objective: to subtract the particle located at peakIdx from the image
% 'cutoutSect'.

% note: 'peakNum' is just to identify which peak is being removed so that
% if the following is executed as a task in parallel it can be easily 
% identified  
peakIdx = sub2ind(size(cutoutSect),peak(1),peak(2),peak(3));

[voronoiVol, Rmask] = peakPlacement(peakIdx,size(cutoutSect),D,w);
cutoutCalcImg = getCalcImg(Rmask,K,D,w,offset,voronoiVol);

sphericalRegion = find(voronoiVol(:) == 1);

diffImg = cutoutSect;
diffImg(sphericalRegion) = abs(cutoutSect(sphericalRegion)-cutoutCalcImg(sphericalRegion));
error2 = diffImg(sphericalRegion).^2;
Rsquared = sum(error2(:));
peakNum = np;
end