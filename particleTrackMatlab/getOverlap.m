function [radiusBottom, heightSphericalCap1, x] = getOverlap(center1,center2,radius1,radius2)
% Objective: 
%   sphere1 and sphere2 intersect at a circle and overlap.  Find the radius
%   of this circle and the height/amt of the overlap region.
% 
% Visual in 2-dimensions:
% 
% case1:
% 
%      (x,y)                            ___
%           / .                          A
%          /    . r2                     |
%       r1/       .                      |y
%        /          .                    |
%      1/_____________.2                _V_
% (0,0) \            . (d,0)
%        \         .
%         \      .
%          \   .
%           \.
%       |<----------->| 
%       |      d      |
%       |             |
%       |<------->|   | 
%       |   r1    |   |
%       |         |   |  
%       |<-->|    |   | 
%       | x  |<-->|   |
%       |     r1-x    |
%       |             |

% y = radius of the circle where sphere1 and sphere2 intersect
% volSphericalCap = ((pi*heightofCap)/6)*(3*y^2 + heightofCap^2)
% 
% volSphericalCap1 = (pi*(r1-x)/6)*(3*(y^2) + (r1-x)^2)
%                  = volume of particle INSIDE overlap region
% 
% volPartice = region OUTSIDE the overlap region. 
% Therefore, volParticle1 = volSphere1 - volSphericalCap1
% 
% ---------OR (case2, x is (-))-------------
% 
%      (x,y)                                    ___
%           \  .                                 A
%            \     . r2                          |
%          r1 \       .                          |y
%              \         .                       |
%          (0,0)\_____________. (d,0)           _V_
%               /
%              / 
%             / 
%            /
%           /
%          
%     |<-------->|            |
%     |    r1    |            |
%     |          |<---------->|
%     |     |<-->|      d     | 
%     |     | x  |     
%     |<--->|
%       r1-x


% volSphericalCap1 = (pi*(r1-x)/6)*(3*(y^2) + (r1-x)^2) 
% ...SAME as CASE1 BUT...
% volSphericalCap1 = volume of particle OUTSIDE overlap region 
% 
% volPartice = region OUTSIDE the overlap region. 
% Therefore, volParticle1 = volSphericalCap1

% IMPORTANT: Don't bother computing volParticle2 because it depends on
% whether (x,y) is to the right (case1) or to the left (case2) of (0,0). It
% ALSO depends on whether r2 < or > d (i.e., if the center of particle2 is
% inside or outside of particle1)

% --------------------------------x--------------------------------
% In reality, the centers of particle 1 and 2 are in 3-Dimensions...
% i.e., 2 spheres intersect at a circle
% (x,y,z) = a point on the circle where the 2 spheres intersect
% 
%      (x,y,z)                          ___
%           / .                          A
%          /    .                        |
%         /       .                      |y
%        /          .                    |
%      1/_____________.2(x2,y2,z2)      _V_
% (x1,y1,z1)
% 
%       |<-->| 
%       |  x    
%       |<------------>|
%              d  

% d = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
d = sqrt((center1(1)-center2(1))^2 + (center1(2)-center2(2))^2 + (center1(3)-center2(3))^2);

% ...to change from 3Dimensions to 2Dimensions...
% 
% Rotate the coord system making sphere 1 the origin and the distance
% between spheres 1 and 2 the NEW x-direction

% Now your only concern is where the great circles of sphere 1 and 2
% intersect at a point (x,y) that is r1 away from the center of sphere 1,
% r2 away from the center of sphere 2, and orthogonal to the NEW
% x-direction between the 2 centers  


% ---------(case1)-------------
% 
%      (x,y,0)                          ___
%           / .                          A
%          / :  . r2                     |
%       r1/  :    .                      |y
%        /  _:      .                    |
%       /__|_:________. (d,0,0)         _V_
% (0,0,0)
% 
%       |<-->|<------>| 
%         x     d - x


% ---------OR (case2)-------------

%    (x,y,0)                                    ___
%           \  .                                 A
%            \     . r2                          |
%          r1 \       .                          |y
%              \         .                       |
%               \_____________. (d,0,0)         _V_
%          (0,0,0)
% 
%           |<-->|<---------->| 
%             x        d

% Eqn of a great circle for sphere 1:
%   (x - 0)^2 + (y-0)^2 = r1^2
%             x^2 + y^2 = r1^2
% 
% Eqn of a great circle for sphere 2:
%   (x - d)^2 + (y-0)^2 = r2^2
%         (x-d)^2 + y^2 = r2^2
% 
% Solve for (x,y)...
%             x^2 + y^2 = r1^2
%         (x-d)^2 + y^2 = r2^2
% 
%          x^2 - (x-d)^2 = r1^2 - r2^2 
%  x^2-(x^2 - 2dx + d^2) = r1^2 - r2^2
%              2dx - d^2 = r1^2 - r2^2
%                    2dx = r1^2 - r2^2 + d^2
%                      x = (r1^2 - r2^2 + d^2)/2d

x = ((radius1^2) - (radius2^2) + (d^2))/(2*d); %x-coord for intersection of spheres 1 and 2. For case1 x is (+); for case2 x is (-)
y = sqrt((radius1^2) - (x^2)); %Note: y is BOTH (-) and (+) for the top and bottom of the circle (i.e., (x,y) and (x,-y))

heightSphericalCap1 = radius1 - abs(x); % r1-x = height of sphericalCap1
radiusBottom = y; % y = radius of the circular bottom of the spherical cap
end