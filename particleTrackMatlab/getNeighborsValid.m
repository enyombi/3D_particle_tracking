function validNeighbors = getNeighborsValid(particlePos,minSepDist)
% objective: to return the positions of particles that do NOT overlap and
% are separated by a distance 'minSepDist'

% particlePos = linear index identifying particle location (row,col,slice)
% dim = size(imageData);
% minSepDist = minimum separation distance between particles for them to 
%              be considered touching

addpath('~/poincareProgs/particleTrackMatlab')

boolPks = getNeighbors(particlePos,minSepDist);

for np = 1:size(boolPks,1)
    %find all peaks that are too close to a given peak 
    %then change the ALL the values of in element for 'peaks_too_close' to 
    %an arbitrary value of 5 so that when the for-loop reaches a column in
    %bool_neighbors that only has 5's it does nothing.
    peaksTooClose = find(boolPks(:,np)==1); 
    boolPks(:,peaksTooClose) = 5;

end

diagonalBool = diag(boolPks);
validPks = find(diagonalBool == 0);
validNeighbors = particlePos(validPks,:); %linear index identifying the positions of particle that don't overlap 
end