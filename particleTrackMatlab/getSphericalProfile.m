function avgProfile = getSphericalProfile(center, D, w, imageData)
% Author: Eru K.
%modified on 1-Oct-2014
%   realized that particleRadius = D/2 (NOT D/2 + 2*w)
%   therefore, 'w' was deleted as input parameter

% Objective: sample the linear profile across a particle at different
% points

% D = particle Diameter
% center = center of particle
% imageParticle = 3Dimage that particle 



% Sphere is circumscribed by a box 
% 
%     slices(+)
%    /
%   /
%   ----->cols(+)
%   |
%   |
%   rows(+)
% 
%    5______________4
%    /:     /'     /|
%   /_:__ 9/_'____/ | 
% 1/:_:__ /:_8___/| |
% |11_: _| c_'_ _||12
% |/:/_ _|/:/____|/ |
% | :   7| /10   ||/ 2
% |/____ |/______|/ 
% 3              6

particleRadius = ((D/2)+(2*w)); %this is the apparent particle radius; modified 1-Oct-2014 
% particleRadius = D/2; %modified 1-Oct-2014, 'w' is a part of D but D is not effected by 'w'

% vertices of diagonal lines through ptC
pts(1,:) = center - particleRadius;
pts(2,:) = center + particleRadius;

pts(3,:) = [center(1)+particleRadius, center(2)-particleRadius, center(3)-particleRadius];
pts(4,:) = [center(1)-particleRadius, center(2)+particleRadius, center(3)+particleRadius];

pts(5,:) = [center(1)-particleRadius, center(2)-particleRadius, center(3)+particleRadius];
pts(6,:) = [center(1)+particleRadius, center(2)+particleRadius, center(3)-particleRadius];

% vertices of vertical and horizontal lines through ptC
pts(7,:) = [center(1), center(2), center(3)-particleRadius];
pts(8,:) = [center(1), center(2), center(3)+particleRadius];

pts(9,:) = [center(1)-particleRadius, center(2), center(3)];
pts(10,:) = [center(1)+particleRadius, center(2), center(3)];

pts(11,:) = [center(1), center(2)-particleRadius, center(3)];
pts(12,:) = [center(1), center(2)+particleRadius, center(3)];

pts = single(pts);
% now make certain that the rectangular box around the particle does not
% exceed the bounds of the 3D image...
for np = 1:length(pts) %particle #
    for dim = 1:3 %each particle has 3 dimensions/coords
        if(pts(np,dim) <= 0)
            pts(np,dim) = 1;
        end

        if(pts(np,dim) > size(imageData,dim))
            pts(np,dim) = size(imageData,dim);
        end
    end
end

% Linear Algebra:
%   ptsP(:,1) = t*A + center(1)
%   ptsP(:,2) = t*B + center(2)
%   ptsP(:,3) = t*C + center(3)
% 
% (A,B,C) = the unit vector for a 3D line from 'center' to pts(np,:) that
% passes through ptsP(np,:)
% 
%   5       3
%     \p5  /
% 1 .  \  /p3
%   p1. c
%     /  \.p2  
%    /p4  \ .
%  4     p6\  .2
%           6     
% 
%   i.e., A, B, C have values such that their distance/magnitude/norm is
%   unity: 
%       A^2 + B^2 + C^2 = 1
% 
%   suppose np = 3,  
%       (ptsP(3,1)-center(1))^2 + (ptsP(3,2)-center(2))^2 + (ptsP(3,3)-center(3))^2 = distance^2
% 
%   then the unit vector from 'center' to 'pts(3,:)' is: 
%       A = (ptsP(:,1)-center(1))/distance
%       B = (ptsP(:,2)-center(2))/distance
%       C = (ptsP(:,3)-center(3))/distance
% 
%   now solve for scalar multiplier t needed for unit vector to go from
%   'center' to ptsP(3,:) i.e., p3
%       constraint: (ptsP(3,1)-center(1))^2 + (ptsP(3,2)-center(2))^2 + (ptsP(3,3)-center(3))^2 = particleRadius^2
% 
%   the equation of the line from 'center' to p3 is given as:
%       t*<A,B,C> + (center(1),center(2),center(3)) = (ptsP(3,1),ptsP(3,2),ptsP(3,3)), vectorize form
% 
%       or in linear algebra form:
%           ptsP(3,1) = t*A + center(1)  ==> t*A = ptsP(3,1) - center(1)
%           ptsP(3,2) = t*B + center(2)  ==> t*B = ptsP(3,2) - center(2)
%           ptsP(3,3) = t*C + center(3)  ==> t*C = ptsP(3,3) - center(3)
% 
%   Substitute, so constraint becomes:
%       (t*A)^2 + (t*B)^2 + (t*C)^2 = particleRadius^2 
%       (t^2)*(A^2 + B^2 + C^2) = particleRadius^2
% 
%   recall, <A,B,C> is a unit vector so (A^2 + B^2 + C^2) = 1
%       t = particleRadius
% 
%   hence,
%       ptsP(3,1) = particleRadius*A + center(1)  
%       ptsP(3,2) = particleRadius*B + center(2)  
%       ptsP(3,3) = particleRadius*C + center(3)  
% 
%   check:
%       constraint: (ptsP(3,1)-center(1))^2 + (ptsP(3,2)-center(2))^2 + (ptsP(3,3)-center(3))^2 = particleRadius^2
%                      (particleRadius*A)^2 + (particleRadius*B)^2    + (particleRadius*C)^2 = particleRadius
%                       particleRadius*(A^2 +          B^2            +                 C^2) = particleRadius
%                                                                             particleRadius = particleRadius                       
% 
% solve for scalar multiplier(t)...
%       (A^2 + B^2 + C^2)*(t^2) = particleRadius^2
%                            t = sqrt(particleRadius^2/(A^2 + B^2 + C^2)
% 
% In summary, you only need to compute <A,B,C> to be able to calculate 
% ptsP(np,:)

ptsP = zeros(size(pts(1:6,:)),'single'); % ptsP(1,:) = pt X; ptsP(2,:) = pt X1; pts(5,:) = y2; pts(6,:) = pt Y
for np = 1:length(ptsP) %pt #
    dist = sqrt((pts(np,1) - center(1))^2 + (pts(np,2) - center(2))^2 + (pts(np,3) - center(3))^2);
    
    A = (pts(np,1) - center(1))/dist; %row/y-directional component of unit vector
    B = (pts(np,2) - center(2))/dist; %column/x-directional component of unit vector
    C = (pts(np,3) - center(3))/dist; %slice/z-directional component of unit vector

    if(dist == 0) %occurs when pts(np,:) == 'center'
        A = 0;
        B = 0;
        C = 0;
    end
    
    ptsP(np,1) = particleRadius.*A + center(1);
    ptsP(np,2) = particleRadius.*B + center(2);
    ptsP(np,3) = particleRadius.*C + center(3);
end

% now make certain that the locations of all ptsP do not
% exceed the bounds of the 3D image...
for np = 1:length(ptsP) %particle #
    for dim = 1:3 %each particle has 3 dimensions/coords
        if(ptsP(np,dim) <= 0)
            ptsP(np,dim) = 1;
        end

        if(ptsP(np,dim) > size(imageData,dim))
            ptsP(np,dim) = size(imageData,dim);
        end
    end
end

% Diagonal profiles:
profileN{1} = getLineProfile(ptsP(1,:),ptsP(2,:),imageData); %profile from P1 to P2
profileN{2} = getLineProfile(ptsP(3,:),ptsP(4,:),imageData); %profile from P3 to P4
profileN{3} = getLineProfile(ptsP(5,:),ptsP(6,:),imageData); %profile from P5 to P6

% Vertical/Horizontal profiles:
profileN{4} = getLineProfile(pts(7,:),pts(8,:),imageData);
profileN{5} = getLineProfile(pts(9,:),pts(10,:),imageData);
profileN{6} = getLineProfile(pts(11,:),pts(12,:),imageData);



% make each profile have the same center and be the same length
% (note: if the  particle is symmetric as it should be, then the number of
% elements in each profile should be the same. The following is just a
% precautionary measure in the event that resolution in x-, y-, z- is not
% the same)
for np = 1:(length(profileN)-1)
    numElements = max(length(profileN{np}),length(profileN{np+1}));
end

profile = zeros([length(profileN), numElements],'single');
for np = 1:length(profileN)
% %     zero padding is NOT good because it introduces trailing zeros that
% %     thwart comparision of shorter and longer profiles. Instead use
% %     linear interpolation to make short profiles that same length while
% % 	keeping the same correlation/function between elements 
%     padsize = floor((numElements - length(profileN{np}))/2);
%     profile(np,:) = padarray(profileN{np}',[0,padsize],0); %padarray(array,[# of rows to pad, # of cols to pad],# to pad with) 
    profile(np,:) = interp1(1:length(profileN{np}),profileN{np},linspace(1,length(profileN{np}),numElements));
end

% lengthening the profiles by linear interpolation should have positioned 
% each profile along the same center so that nearly all profiles have the 
% same values 
avgProfile = mean(profile,1);
end

% % Testing getSphericalProfile.m
% D = 2;
% center = [2 2 2];
%  
% imageData(:,:,1) = magic(3);
% imageData(:,:,2) = magic(3)+1;
% imageData(:,:,3) = magic(3)-1;

% -------------
% Test using getSphericalProfile.m to get the spherical profile and use
% nlinfit() with ssf.m on this to get paramters K, D, and w.
% K = 0.9;
% D = 120;
% w = 12;
% 
% ss=2*fix(D/2+4*w/2)-1;         % size of ideal particle image
% os=(ss-1)/2;                   % (size-1)/2 of ideal particle image
% center = [os+1 os+1 os+1]; %center of ideal particle image
% [voronoiVol Rmask] = peakPlacement(center,[ss ss ss],D,w);
% imageData = getCalcImg(Rmask,K,D,w,voronoiVol);
% 
% avgProfile = getSphericalProfile(center,D,w,imageData);

% b = ceil((length(avgProfile)-1)/2);
% r = -b:b;
% param = nlinfit(r,avgProfile,'ssf',[1 110 10]); 
% %values returned: 
% % K = 0.9000 (acutal 0.9)
% % D = 119.8995 (actual 120)
% % w = 12.0069 (actual 12)






% ------------------ (math bkgrd)--------------------

% Eqn of a 1D vector: x = t*a+x0 
%     where,
%       a = directional component
%       t = scalar multiplier
%       x0 = starting value 
% 
% i.e., x = [x0  _  _  _  _ ...] 
%                -->a*t
% 
% A line traditional has 2 directions: x,y.  Therefore,
%   x = t*A + x0
%   y = t*B + y0
% 
% (x,y) = t*<A,B> + (x0,y0)
% 
% In 3Dimensions, a line is:
%   (x,y,z) = t*<A,B,C> + (x0,y0,z0)
% 
% Find pts x, x1, y, and y1 by calculating their intersection with the
% sphere: x^2 + y^2 + z^2 = (D/2)^2

% Linear Algebra:
%   x = t*A + x0
%   y = t*B + y0
%   z = t*C + z0
% 
%   constraint: (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = (D/2)^2
% 
% let x = row, y = col, and z = slice 
% therefore, (x0,y0,z0) = [center(1) center(2) center(3)]



% ------------------ (End: math bkgrd)--------------------
