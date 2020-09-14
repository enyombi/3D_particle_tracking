function [sepDist, touchingBool] = getSepDist(center1,distance,center2,D2)
% Objective: To determine if a sphere at center1 with radius 'distance'
% intersects a sphere with radius (D2/2) at center2 

% center1, center2 = position of particles1,2 in 3D image [row, col, slice]
% distance = distance from 'center1'
% D2 = size/diameter of particle2 

% step1: compute the vector/(direction of line) from center1 to center2
% step2: find the location of a point M (xm,ym,zm) located at a distance
%        'dist' from 'center1'
% step3: Determine if pt. M (xm,ym,zm) is close enough to 'center2' for it
% to be touching [i.e., determine if (xm,ym,zm) is located a distance <=
% D2/2]

% step1:
diff = center2 - center1;
magnitude = sqrt(sum(diff.^2)); %distance btwn center1 and center2
vector = diff./magnitude; %unit-vector from center1 to center2

% A = diff(1)/magnitude; %row-direction (y)
% B = diff(2)/magnitude; %column-direction (x)
% C = diff(3)/magnitude; %slice-direction (z)

% step2: 
posM = (vector.*distance) + center1; %position of pt. M

% step3: 
dist2M = sqrt(sum((center2 - posM).^2)); %distance btwn center2 to pt. M
edge2 = (-1*vector)*(D2/2) + center2; %position of a point on the line from center2 to center1 (i.e., -1*vector) that is on the edge of the sphere at center2
sepDist = sqrt(sum((posM-edge2).^2)); %distance between pt. M and the edge of the sphere at center2 along vector from center1 to center2

% CRUCIAL NOTE:
% 'sepDist' is a scalar that is the distance btwn pt. M and the edge of te
% sphere at center2
% 
% if the sphere at center2 overlaps the sphere of radius 'distance'
% surrounding center1 (i.e., touchingBool = 1) then 'sepDist' == amount of
% overlap!  

touchingBool = 0;
if(dist2M <= (D2/2)) %(D2/2) = radius of sphere at center2
    touchingBool = 1;
end
end