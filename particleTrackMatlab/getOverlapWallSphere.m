function [radiusBottom, heightSphericalCap, x] = getOverlapWallSphere(center,radius,dimC,dim,boolLowerBound)
% Objective: 
%   a sphere intersects a wall/plane at a circle

% ----------------------------------x
% INPUT vars:
% center = center of sphere
% radius = radius of sphere
% dimC(dim) = boundary closest to center(dim).  All other values
%             in dimC are irrelevant
% boolLowerBound = 1 if dimC(dim) is one of the lower bounds (See e.g2)

% e.g1, 
%   if sphere's center = [3.5,2,2] & radius = 1 in a 4x4x4 box
% 
%   and you want to find the part of the sphere sticking out of the right
%   face of the box (i.e., dim = 1)
% 
%   then, dimC = [4,nan,nan]

% e.g2, 
%   if sphere's center = [2.5,2,2] & radius = 1 in a 4x4x4 box
% 
%   and you want to find the part of the sphere sticking out of the left
%   face of the box (i.e., dim = 1)
% 
%   then, dimC = [0,nan,nan]

% x < 0 even though center is inside the box.
% Therefore, let boolLowerBound ==1 in order to reverse the sign of 'x'
% ----------------------------------x

% OUTPUT vars:
%   radiusBottom = radius of the circular intersection
%   heightSphericalCap = height of the sphere that crosses over the wall
%       Note: if heightSphericalCap < 0, then sphere doesn't cross over the
%       wall
%   x = dimC(dim) - center(dim), where dim = 1, 2, or 3
%    Note (this only applies when dimC(dim) >> 0):
%     if x > 0, then center is INSIDE the container (case1)
%     if x = 0, then center is on the wall (case1)
%     if x < 0, then center is OUTSIDE of the container (case2)

%   abs(x) = distance between the center of the sphere and the wall in a
%            single dimension 
% 

% for dim = 1, find where the sphere interesects the wall in the
% row-dimension: 

% case1:
% 
%                   (dimC(1),center(2)+radiusBottom,center(3))          ___
%                                    /.                                  A
%                                   / .                                  |
%                          radius  /  .                      radiusBottom|
%                                 /   .                                  |
% (center(1),center(2),center(3))/____._____(center(1)+radius,center(2))_V_
%                                \    .(dimC(1),center(2),center(3))
%                                 \   .
%                                  \  .
%                                   \ .
%                                    \.
%                            (dimC(1),center(2)-radiusBottom,center(3))
%                               |     |     |
%                               |<--------->|
%                               |   radius  |
%                               |<--->|     |
%                                 |x| |<--->|
%                                     |  h   
%                                     |
%                          |---------------------->(+)x-direction 
%                                     |                              
%                                     |  
%                                     |--->(OUTSIDE container)
%               (INSIDE container)<---|
% 
% x = dimC(1) - center(1)
% h = radius - |x|; %heightofCap
% 
% solve for radiusBottom using Pythagorean Theorem:
%                     radius^2 = |x|^2 + radiusBottom^2
%       (radius^2 - |x|^2) = radiusBottom^2
%   sqrt(radius^2 - |x|^2) = radiusBottom

x = dimC(dim) - center(dim);

heightSphericalCap = radius - abs(x); %(center(dim)+radius) - dimC(dim)

radiusBottom = sqrt((radius^2) - (abs(x)^2)); %vertical distance above (dimC(1),center(2),center(2)) where the sphere intersects the container which is ALSO the radius of the circle where the sphere intersects the wall


% x > 0. So center is inside of the container.  Only sliver whose volume is
% 'volSphericalCap' of the particle is outside of the container.
% 
% volSphericalCap = ((pi*heightofCap)/6)*(3*y^2 + heightofCap^2)
%                 = ((pi*h)/6)*(3*radiusBottom^2 + h^2)
% 
% Therefore, volParticle = volSphere - volSphericalCap
% 

% % do this outside this function:
% yPos = center(2) + radiusBottom;
% yNeg = center(2) - radiusBottom;
% 
% ptsInteresect(1,:) = [x,yPos];
% ptsInteresect(1,:) = [x,yNeg];




% case2 (center outside the container walls):
% 
%                   (dimC(1),y,center(3))                           ___
%                                .\                                  A
%                                . \                                 |
%                                .  \radius              radiusBottom|
%   (dimC(1),center(2),center(3)).   \                               |
% (center(1)-radius,center(2))___.____\                             _V__ 
%                                      (center(1),center(2),center(3))
%                                |<--->| 
%                                | |x| |
%                            |<------->|
%                            |  radius
%                                |
%                            |<->|
%                             radius-|x|
%                                |   
%                                |------>(+)x-direction 
%                                |
%                                |--->(OUTSIDE container)
%        (INSIDE container)<-----|
% 
% x = dimC(1)-center(1); 
%   dimC(1) < center(1) so, x < 0. This means most of the particle is
%   outside the container. Only the sphericalCap is inside
% 
% h = radius - |x| 
% 
% volSphericalCap = ((pi*heightofCap)/6)*(3*y^2 + heightofCap^2)
% Therefore, volParticle = volSphericalCap

if(boolLowerBound == 1) %boolLowerBound = 1 or boolUpperBound = 0
    x = -1*x;
end
end