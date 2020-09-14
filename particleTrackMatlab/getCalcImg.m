function ipi = getCalcImg(r,K,D,w,offset,voronoiVol)
% objective: To draw the image of particle ONLY in the spherical region of
% the 3D cube occupied by the particle.

% adopted from ipf_eru.m

%'voronoiVol' is an image where the pixels for each particle are assigned
%that particle's number (i.e, all the pixels for particle 1 are assigned 1,
%all the pixels for particle 2 are 2, etc....so max(voroni_vol(:)) is the
%number of the last particle in the image.)

addpath('~/poincareProgs/particleTrackMatlab')

calcImg = min(offset(:))*ones(size(voronoiVol),'single'); %set all pixel values in background to the lowest offset value found for all particles

for np = 1:max(voronoiVol(:))
    sphericalRegion = find(voronoiVol == np);
    calcImg(sphericalRegion) = ssf([K(np), D(np), w(np), offset(np)],r(sphericalRegion)); %the for-loop draws one particle at a time in SPHERICAL regions of the image 

%     disp(sprintf('\nIn getCalcImg.m...finished drawing particle at peak %d of %d\n',np,max(voronoiVol(:))))
    clear sphericalRegion
end

ipi = calcImg;
end