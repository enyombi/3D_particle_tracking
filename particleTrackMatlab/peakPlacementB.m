function [voronoiVol, r] = peakPlacementB(centers,dim,D,w,dp)
%Author: Eru Kyeyune-Nyombi
% filename: peakPlacementB

% adopted from and exactly the same as peakPlacement.m except recieves the
% subscripted (i.e., row,column,slice) version of centers and is therefore
% capable of position particle at sub-pixel locations

% Objective: To return the SPHERICAL region of radial values used to
% calculate the ideal mask particle image

% adopted from peak_placement.m. This is a trimmed down version of
% peak_placement.m that uses the actual particle size, i.e. ((D/2)+(2*w)),
% to define spherical region occupied by the particle

% centers = particle positions 
% dim = dimensions of the 3D image dim(1) = # of rows, dim(2) = # of cols,
% etc.
% D = particle diameter
% w = transition region at edges of the particle
% dp = displacement (used for locating the center at sub-pixel resolution
% by adding/subtracting <1 pixel to row, col, slice)
%   r = sqrt((row - center(1) + dp(1)).^2 + (col - center(2) + dp(2)).^2 + (slice - center(3) + dp(3)).^2);}

dim = single(dim);
[row, col, slice] = ndgrid(1:dim(1),1:dim(2),1:dim(3));

if ~exist('dp','var')
    row = single(row);
    col = single(col);
    slice = single(slice);
end

if exist('dp','var')
    row = single(row)+dp(1);
    col = single(col)+dp(2);
    slice = single(slice)+dp(3);
end

% bkgrdVal = max(dim);%arbitrary value used so that when ipf([D w], r) is called the particle is distinct from the bkgrd (note: bkgrdVal should be outside the range of values used in constructing the mask particle image) 
r = max(dim)*ones(dim,'single'); %radial values for the calculated image
voronoiVol = zeros(dim,'single'); %eqivalent to 'over' in Shattuck's 'pgrid.m'

% -------------------(start: modification on 22Aug2014)--------------------
% 'voronoiVol' identifies the space occupied by each particle in the image.
% These spaces do NOT overlap!  And the edges of particles defined by 'w'
% have very small pixel values so that if the space alotted to one particle
% with a very large w overlaps a particle.  And this particle with the
% large w is deleted first. Then a large portion of the overlapped
% paritlcle is likely to remain in the image. 

% Therefore, to ensure a particle is entirely deleted, prioritize the
% assignment of spaces to particles in 'voronoiVol' with the particles
% having the lowest ratio w/D being drawn first so that they are completely
% deleted from the image.

% correction: you want the particle with large w to be drawn in
% 'voronoiVol' first.  This way the particle it overlaps will be deleted
% last...Think overlapping spheres in 'voronoiVol'. The only way the
% particle 578 will be COMPLETELY deleted is if it is drawn AFTER particle
% 901.

ratio = w./D;
[naught,order] = sort(ratio,'descend'); %sort in ascending order
% -------------------(end: modification on 22Aug2014)--------------------

for np = 1:size(centers,1)
%     rowTemp = row - centers(np,1); %y-displacement from center(np, )
%     colTemp = col - centers(np,2); %x-displacement from center(np, )
%     sliceTemp = slice - centers(np,3); %z-displacment from center(np, )

    %Eqn of Sphere:
    %     (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = Radius^2
    %     where center = (x0,y0,z0)
    rTemp = sqrt((row - centers(order(np),1)).^2 + (col - centers(order(np),2)).^2 + (slice - centers(order(np),3)).^2); 
             
    RadiusApparent = ((D(order(np))/2)+(2*w(order(np))));
    sphericalRegion = find(rTemp(:) <= RadiusApparent); 
    
    r(sphericalRegion) = rTemp(sphericalRegion);
    
    voronoiVol(sphericalRegion) = order(np);
%     disp(sprintf('In peakPlacement.m...finished locating peak %d of %d\n',order(np),size(centers,1)))
end
end